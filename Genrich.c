/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.1
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <zlib.h>
#include <omp.h>
#include "Genrich.h"

/* void printVersion()
 * Print version and copyright.
 */
void printVersion(void) {
  fprintf(stderr, "Genrich, version %s\n", VERSION);
  fprintf(stderr, "Copyright (C) 2018 John M. Gaspar (jsh58@wildcats.unh.edu)\n");
  exit(-1);
}

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./Genrich  -%c <file>  -%c <file>", INFILE, OUTFILE);
  fprintf(stderr, "  [optional arguments]\n");
  fprintf(stderr, "Required arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input SAM/BAM file(s) for treatment sample(s)\n", INFILE);
  fprintf(stderr, "  -%c  <file>       Output peak file\n", OUTFILE);
  fprintf(stderr, "Optional arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input SAM/BAM file(s) for control sample(s)\n", CTRLFILE);
  fprintf(stderr, "  -%c  <file>       Output bedgraph file for pileups and p/q values\n", LOGFILE);
  fprintf(stderr, "Filtering options:\n");
  fprintf(stderr, "  -%c  <arg>        Comma-separated list of chromosomes to ignore\n", XCHROM);
  fprintf(stderr, "  -%c  <int>        Minimum MAPQ to keep an alignment (def. 0)\n", MINMAPQ);
  fprintf(stderr, "Options for unpaired alignments:\n");
  fprintf(stderr, "  -%c               Keep unpaired alignments (def. false)\n", SINGLEOPT);
  fprintf(stderr, "  -%c  <int>        Keep unpaired alignments, with fragment length\n", EXTENDOPT);
  fprintf(stderr, "                     increased to specified value\n");
  fprintf(stderr, "  -%c               Keep unpaired alignments, with fragment length\n", AVGEXTOPT);
  fprintf(stderr, "                     increased to average value of paired alignments\n");
  fprintf(stderr, "Options for peak calling:\n");
  fprintf(stderr, "  -%c  <float>      Maximum q-value (FDR-adjusted p-value; def. %.2f)\n", QVALUE, DEFQVAL);
  fprintf(stderr, "  -%c  <float>      Maximum p-value (mutually exclusive with -%c)\n", PVALUE, QVALUE);
  fprintf(stderr, "  -%c  <int>        Minimum length of a peak, in bp (def. %d)\n", MINLEN, DEFMINLEN);
  fprintf(stderr, "  -%c  <int>        Maximum distance between significant sites to\n", MAXGAP);
  fprintf(stderr, "                     cluster them together, in bp (def. %d)\n", DEFMAXGAP);
  fprintf(stderr, "I/O options:\n");
  fprintf(stderr, "  -%c               Option to gzip-compress output(s)\n", GZOPT);
  fprintf(stderr, "  -%c               Option to print status updates/counts to stderr\n", VERBOSE);
/*
  fprintf(stderr, "  -%c  <int>        Number of threads to use (def. %d)\n", THREADS, DEFTHR);
*/
  exit(-1);
}

/*** Utilites ***/

/* int error()
 * Prints an error message.
 */
int error(char* msg, enum errCode err) {
  fprintf(stderr, "Error! %s%s\n", msg, errMsg[err]);
  return -1;
}

/* void* memalloc()
 * Allocates a heap block.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* void* memrealloc()
 * Changes the size of a heap block.
 */
void* memrealloc(void* ptr, int size) {
  void* ans = realloc(ptr, size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

/* char* getLine()
 * Reads the next line from a file.
 */
char* getLine(char* line, int size, File in, bool gz) {
  if (gz)
    return gzgets(in.gzf, line, size);
  else
    return fgets(line, size, in.f);
}

/*** Quicksort ***/
// adapted from https://www.geeksforgeeks.org/quick-sort/

/* void swapFloat(): Swap two float values (pileup->cov)
 * void swapInt():   Swap two int values (pileup->end)
 * int partition():  Place last elt into correct spot
 * void quickSort(): Control quickSort process recursively
 */
void swapFloat(float* a, float* b) {
  float t = *a;
  *a = *b;
  *b = t;
}
void swapInt(unsigned int* a, unsigned int* b) {
  unsigned int t = *a;
  *a = *b;
  *b = t;
}
int partition(Pileup* p, int low, int high) {
  float pivot = p->cov[high];  // pivot value: last elt
  int idx = low - 1;

  for (int j = low; j < high; j++) {
    if (p->cov[j] < pivot) {
      idx++;
      swapFloat(&p->cov[idx], &p->cov[j]);
      swapInt(&p->end[idx], &p->end[j]);  // swap int values too
    }
  }
  idx++;
  swapFloat(&p->cov[idx], &p->cov[high]);
  swapInt(&p->end[idx], &p->end[high]);
  return idx;
}
void quickSort(Pileup* p, int low, int high) {
  if (low < high) {
    int idx = partition(p, low, high);
    quickSort(p, low, idx - 1);
    quickSort(p, idx + 1, high);
  }
}

/*** Calculate q-values ***/

/* float lookup()
 * Return the pre-computed q-value for a given p-value,
 *   using parallel arrays (p->cov and qval).
 */
float lookup(Pileup* p, int low, int high, float* qval,
    float pval) {
  if (low == high)
    return qval[low];
  int idx = (low + high) / 2;
  if (p->cov[idx] == pval)
    return qval[idx];
  if (p->cov[idx] > pval)
    return lookup(p, low, idx - 1, qval, pval);
  return lookup(p, idx + 1, high, qval, pval);
}

/* void saveQval()
 * Calculate and save q-values, given the pre-compiled
 *   list of p-values (Pileup* p).
 */
void saveQval(Chrom* chrom, int chromLen,
    unsigned long genomeLen, Pileup* p, int pLen) {

  // sort pileup by p-values
  quickSort(p, 0, pLen - 1);

  // calculate q-values for each p-value: -log(q) = -log(p*N/k)
  unsigned int k = 1;  // 1 + number of bases with higher -log(p)
  float logN = -1 * log10f(genomeLen);
  float* qval = (float*) memalloc((pLen + 1) * sizeof(float));
  qval[pLen] = FLT_MAX;
  for (int i = pLen - 1; i > -1; i--) {
    // ensure monotonicity
    qval[i] = MAX( MIN( p->cov[i] + logN + log10f(k),
      qval[i + 1]), 0.0f);
    k += p->end[i];
  }

  // save pileups of q-values for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;
    for (int j = 0; j < chr->pvalLen; j++)
      chr->qval->cov[j] = lookup(p, 0, pLen, qval,
        chr->pval->cov[j]);
  }

  // free memory
  free(qval);
  free(p->end);
  free(p->cov);
  free(p);
}

/*** Calculate p-values ***/

/* void saveConst()
 * Save a given value as pileup for a full chromosome.
 */
void saveConst(Pileup** p, int* size, int len, float val) {
  (*p) = (Pileup*) memalloc(sizeof(Pileup));
  (*p)->end = (unsigned int*) memalloc(sizeof(unsigned int));
  (*p)->end[0] = len;
  (*p)->cov = (float*) memalloc(sizeof(float));
  (*p)->cov[0] = val;
  *size = 1;
}

/* void recordPval()
 * Save length of given p-value into Pileup arrays.
 */
void recordPval(Pileup* p, int* pLen, int* pMem,
    float pval, int length) {
  // look for p-value in array
  for (int i = 0; i < *pLen; i++)
    if (p->cov[i] == pval) {
      // if found, just add length
      p->end[i] += length;
      return;
    }

  // realloc memory if necessary
  if (*pLen + 1 > *pMem) {
    *pMem += 10000;
    p->end = (unsigned int*) memrealloc(p->end,
      *pMem * sizeof(unsigned int));
    p->cov = (float*) memrealloc(p->cov,
      *pMem * sizeof(float));
  }

  // save values to arrays
  p->end[*pLen] = length;
  p->cov[*pLen] = pval;
  (*pLen)++;
}

/* int countIntervals()
 * Count the number of pileup intervals to create
 *   for a composite.
 */
int countIntervals(Chrom* chr) {
  int num = 0;
  int k = 0;
  for (int j = 0; j < chr->treatLen; j++) {
    while (k < chr->ctrlLen
        && chr->ctrl->end[k] < chr->treat->end[j]) {
      num++;
      k++;
    }
    if (chr->ctrl->end[k] == chr->treat->end[j])
      k++;
    num++;
  }
  return num;
}

/* float calcPval()
 * Calculate -log10(p) using an exponential distribution
 *   with parameter beta ('ctrl') and observation 'treat'.
 */
float calcPval(float treatVal, float ctrlVal) {
  if (ctrlVal == 0.0f)
    return treatVal == 0.0f ? 0.0f : FLT_MAX;
  return treatVal / ctrlVal * M_LOG10E;
}

/* Pileup* savePval()
 * Create and save p-values as pileups for each chrom.
 *   Return a pileup of all p-values to be used in
 *   q-value calculations.
 */
Pileup* savePval(Chrom* chrom, int chromLen,
    bool qvalOpt, int* pLen) {

  // create pileup for conversion of p-values to q-values
  Pileup* p = NULL;
  int pMem = 0;
  if (qvalOpt) {
    p = (Pileup*) memalloc(sizeof(Pileup));
    p->end = (unsigned int*) memalloc(10000 * sizeof(unsigned int));
    p->cov = (float*) memalloc(10000 * sizeof(float));
    pMem = 10000;
  }

  // create pileups for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // check pileups for missing values
    if (chr->treat == NULL)
      saveConst(&chr->treat, &chr->treatLen, chr->len, 0.0f);
    if (chr->ctrl == NULL)
      exit(error(chr->name, ERRCTRL));

    // create 'pileup' arrays for p-values
    int num = countIntervals(chr);
    chr->pval = (Pileup*) memalloc(sizeof(Pileup));
    chr->pval->end = (unsigned int*) memalloc(num * sizeof(unsigned int));
    chr->pval->cov = (float*) memalloc(num * sizeof(float));
    if (qvalOpt) {
      // create pileups arrays for q-values (to be populated later)
      chr->qval = (Pileup*) memalloc(sizeof(Pileup));
      chr->qval->end = (unsigned int*) memalloc(num * sizeof(unsigned int));
      chr->qval->cov = (float*) memalloc(num * sizeof(float));
    }
    chr->pvalLen = num;

    // save p-values to arrays
    int start = 0;
    int j = 0, k = 0;
    for (int m = 0; m < num; m++) {
      if (chr->ctrl->end[k] < chr->treat->end[j]) {
        chr->pval->end[m] = chr->ctrl->end[k];
        chr->pval->cov[m] = calcPval(chr->treat->cov[j],
          chr->ctrl->cov[k]);
        k++;
      } else {
        chr->pval->end[m] = chr->treat->end[j];
        chr->pval->cov[m] = calcPval(chr->treat->cov[j],
          chr->ctrl->cov[k]);
        if (chr->ctrl->end[k] == chr->treat->end[j])
          k++;
        j++;
      }
      if (qvalOpt) {
        recordPval(p, pLen, &pMem, chr->pval->cov[m],
          chr->pval->end[m] - start);
        start = chr->pval->end[m];
      }
    }

  }
  return p;
}

/*** Call peaks ***/

/* void printInterval()
 * Print bedgraph(ish) interval, with values of:
 *   pileups (treatment and control), -log(p),
 *   -log(q), and significance ('*') for each.
 */
void printInterval(File out, bool gzOut, char* name,
    int start, int end, float treatVal, float ctrlVal,
    float pval, float qval, bool sig) {
  if (gzOut) {
    gzprintf(out.gzf, "%s\t%d\t%d\t%.5f\t%.5f\t%.5f",
      name, start, end, treatVal, ctrlVal, pval);
    if (qval != -1.0f)
      gzprintf(out.gzf, "\t%.5f", qval);
    gzprintf(out.gzf, "%s\n", sig ? "\t*" : "");
  } else {
    fprintf(out.f, "%s\t%d\t%d\t%.5f\t%.5f\t%.5f",
      name, start, end, treatVal, ctrlVal, pval);
    if (qval != -1.0f)
      fprintf(out.f, "\t%.5f", qval);
    fprintf(out.f, "%s\n", sig ? "\t*" : "");
  }
}

/* void printPeak()
 * Print BED interval for a peak.
 */
void printPeak(File out, bool gzOut, char* name, int start,
    int end, int count, float val, float fe,
    float pval, float qval, int pos) {
  if (gzOut)
    gzprintf(out.gzf, "%s\t%d\t%d\tpeak_%d\t%d\t.\t%f\t%f\t%f\t%d\n",
      name, start, end, count,
      MIN((unsigned int) (val * 10.0f + 0.5f), 1000),
      fe, pval, qval, pos);
  else
    fprintf(out.f, "%s\t%d\t%d\tpeak_%d\t%d\t.\t%f\t%f\t%f\t%d\n",
      name, start, end, count,
      MIN((unsigned int) (val * 10.0f + 0.5f), 1000),
      fe, pval, qval, pos);
}

/* int callPeaks()
 * Call peaks. Produce output on the fly.
 *   Return number of peaks.
 */
int callPeaks(File out, File log, bool logOpt, bool gzOut,
    Chrom* chrom, int chromLen, float pqvalue,
    bool qvalOpt, int minLen, int maxGap) {

  // loop through chroms
  int count = 0;      // count of peaks
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // reset peak variables
    int peakStart = -1, peakEnd = -1; // ends of potential peak
    float summitVal = -1.0f;          // summit p/q value
    int summitPos = 0;                // distance from peakStart to summit
    int summitLen = 0;                // length of summit interval
    float summitPval = -1.0f, summitQval = -1.0f; // summit p- and q-values
    float summitFE = -1.0f;           // summit fold enrichment

    // loop through intervals (defined by chr->pval)
    int start = 0;      // start of interval
    int j = 0, k = 0;   // indexes into chr->treat, chr->ctrl
    for (int m = 0; m < chr->pvalLen; m++) {

      bool sig = false;
      float val = qvalOpt ? chr->qval->cov[m] : chr->pval->cov[m];
      if ( val > pqvalue ) {

        // interval reaches significance
        sig = true;
        if (peakStart == -1)
          peakStart = start;  // start new potential peak
        peakEnd = chr->pval->end[m];  // end of potential peak

        // check if interval is summit for this peak
        if (val > summitVal) {
          summitVal = val;
          summitFE = chr->ctrl->cov[k] ?
            chr->treat->cov[j] / chr->ctrl->cov[k] : FLT_MAX;
          summitPval = chr->pval->cov[m];
          summitQval = qvalOpt ? chr->qval->cov[m] : -1.0f;
          summitPos = (peakEnd + (m ? chr->pval->end[m-1] : 0) ) / 2
            - peakStart;  // midpoint of interval
          summitLen = peakEnd - (m ? chr->pval->end[m-1] : 0);
        } else if (val == summitVal) {
          // update summitPos only if interval is longer
          int len = chr->pval->end[m] - (m ? chr->pval->end[m-1] : 0);
          if (len > summitLen) {
            summitPos = (peakEnd + (m ? chr->pval->end[m-1] : 0) ) / 2
              - peakStart;  // midpoint of interval
            summitLen = len;
          }
        }

      } else {

        // interval does not reach significance
        if (peakStart != -1
            && chr->pval->end[m] - peakEnd > maxGap) {
          // determine if prior peak meets length threshold
          if (peakEnd - peakStart >= minLen) {
            printPeak(out, gzOut, chr->name, peakStart,
              peakEnd, count, summitVal, summitFE,
              summitPval, summitQval, summitPos);
            count++;
          }
          peakStart = -1;     // reset peak start
          summitVal = -1.0f;  // reset peak summit value
          summitLen = 0;      // reset peak summit length
        }
      }

      if (logOpt)
        // print all stats for interval
        printInterval(log, gzOut, chr->name,
          start, chr->pval->end[m],
          chr->treat->cov[j], chr->ctrl->cov[k],
          chr->pval->cov[m],
          qvalOpt ? chr->qval->cov[m] : -1.0f, sig);

      // update chr->treat and chr->ctrl indexes
      if (chr->ctrl->end[k] < chr->treat->end[j])
        k++;
      else {
        if (chr->ctrl->end[k] == chr->treat->end[j])
          k++;
        j++;
      }

      start = chr->pval->end[m];
    }

    // determine if last peak meets length threshold
    if (peakStart != -1 && peakEnd - peakStart >= minLen) {
      printPeak(out, gzOut, chr->name, peakStart,
        peakEnd, count, summitVal, summitFE,
        summitPval, summitQval, summitPos);
      count++;
    }

  }

  return count;
}

/* void findPeaks()
 * Control process of finding peaks:
 *   calculating p- and q-values, calling peaks,
 *   and printing output.
 */
void findPeaks(File out, File log, bool logOpt, bool gzOut,
    unsigned long genomeLen, Chrom* chrom, int chromLen,
    float pqvalue, bool qvalOpt, int minLen, int maxGap,
    bool verbose) {

  if (verbose) {
    fprintf(stderr, "Peak-calling parameters:\n");
    fprintf(stderr, "  Genome length: %ldbp\n", genomeLen);
    fprintf(stderr, "  Significance threshold: -log(%c) > %.3f\n",
      qvalOpt ? 'q' : 'p', pqvalue);
    fprintf(stderr, "  Max. gap between sites: %dbp\n", maxGap);
    fprintf(stderr, "  Min. peak length: %dbp\n", minLen);
  }

  // compute p- and q-values
  int pLen = 0;
  Pileup* p = savePval(chrom, chromLen, qvalOpt, &pLen);
  if (qvalOpt)
    // convert p-values to q-values
    saveQval(chrom, chromLen, genomeLen, p, pLen);

  // identify peaks
  int count = callPeaks(out, log, logOpt, gzOut, chrom,
    chromLen, pqvalue, qvalOpt, minLen, maxGap);
  if (verbose)
    fprintf(stderr, "Peaks identified: %d\n", count);

}

/*** Save treatment/control pileup values ***/

/* float calcLambda()
 * Calculate a background lambda value: sum of fragment
 *   lengths divided by total genome length.
 */
float calcLambda(Chrom* chrom, int chromLen,
    double fragLen, unsigned long* genomeLen) {
  for (int i = 0; i < chromLen; i++)
    if (! (chrom + i)->skip)
      *genomeLen += (chrom + i)->len;
  if (! *genomeLen)
    exit(error("", ERRGEN));
  return fragLen / *genomeLen;
}

/* void savePileupNoCtrl()
 * When no control is available, save the control
 *   pileup as the background lambda value.
 */
void savePileupNoCtrl(Chrom* chrom, int chromLen,
    double fragLen, unsigned long* genomeLen) {
  float lambda = calcLambda(chrom, chromLen, fragLen,
    genomeLen);
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;
    saveConst(&chr->ctrl, &chr->ctrlLen, chr->len, lambda);
  }
}

/* void savePileupCtrl()
 * Save pileup values for control sample(s) from
 *   'diff' arrays and background lambda value.
 */
void savePileupCtrl(Chrom* chrom, int chromLen,
    double fragLen, unsigned long* genomeLen,
    float epsilon) {

  // calculate background lambda value
  float lambda = calcLambda(chrom, chromLen, fragLen,
    genomeLen);

  // create pileup for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // if no read coverage, save constant pileup of 0
    if (chr->diff == NULL) {
      saveConst(&chr->ctrl, &chr->ctrlLen, chr->len, 0.0f);
      continue;
    }

    // determine number of pileup intervals
    int num = 1;
    for (int j = 1; j < chr->len; j++)
      if (fabsf(chr->diff[j]) > epsilon)
        num++;

    // create pileup arrays
    chr->ctrl = (Pileup*) memalloc(sizeof(Pileup));
    chr->ctrl->end = (unsigned int*) memalloc(num * sizeof(unsigned int));
    chr->ctrl->cov = (float*) memalloc(num * sizeof(float));
    chr->ctrlLen = num;

    // save pileup values
    int pos = 0;              // position in pileup arrays
    float val = chr->diff[0]; // pileup value
    int j;
    for (j = 1; j < chr->len; j++)
      if (fabsf(chr->diff[j]) > epsilon) {
        chr->ctrl->end[pos] = j;
        chr->ctrl->cov[pos] = MAX(val, lambda);

        // update pileup value
        val += chr->diff[j];
        // reset val to nearest int if it's close
        if (val < 0.0f)
          val = 0.0f;
        else {
          int valInt = (int) (val + 0.5f);
          if (fabsf(valInt - val) < epsilon)
            val = (float) valInt;
        }

        pos++;
      }

    // save final interval
    chr->ctrl->end[pos] = j;
    chr->ctrl->cov[pos] = MAX(val, lambda);
  }

}

/* double savePileupTreat()
 * Save pileup values for treatment sample(s) from
 *   'diff' arrays.
 *   Return total length of all fragments (weighted).
 */
double savePileupTreat(Chrom* chrom, int chromLen,
    float epsilon) {

  // create pileup for each chrom
  double fragLen = 0.0;  // weighted fragment length
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // if no read coverage, save constant pileup of 0
    if (chr->diff == NULL) {
      saveConst(&chr->treat, &chr->treatLen, chr->len, 0.0f);
      continue;
    }

    // determine number of pileup intervals
    int num = 1;
    for (int j = 1; j < chr->len; j++)
      if (fabsf(chr->diff[j]) > epsilon)
        num++;

    // create pileup arrays
    chr->treat = (Pileup*) memalloc(sizeof(Pileup));
    chr->treat->end = (unsigned int*) memalloc(num * sizeof(unsigned int));
    chr->treat->cov = (float*) memalloc(num * sizeof(float));
    chr->treatLen = num;

    // save pileup values
    int pos = 0;              // position in pileup arrays
    float val = chr->diff[0]; // pileup value
    int start = 0;            // beginning coordinate of interval
    int j;
    for (j = 1; j < chr->len; j++)
      if (fabsf(chr->diff[j]) > epsilon) {
        chr->treat->end[pos] = j;
        chr->treat->cov[pos] = val;

        // keep track of fragment length (weighted by val)
        fragLen += (j - start) * val;
        start = j;

        // update pileup value
        val += chr->diff[j];
        // reset val to nearest int if it's close
        if (val < 0.0f)
          val = 0.0f;
        else {
          int valInt = (int) (val + 0.5f);
          if (fabsf(valInt - val) < epsilon)
            val = (float) valInt;
        }

        pos++;
      }

    // save final interval
    chr->treat->end[pos] = j;
    chr->treat->cov[pos] = val;
    fragLen += (j - start) * val;
  }

  if (fragLen == 0.0)
    exit(error("", ERRTREAT));
  return fragLen;
}

/*** Convert alignments to intervals ***/

/* int saveInterval()
 * Save BED interval for a read/fragment.
 *   Return fragment length.
 */
int saveInterval(Chrom* chrom, int start, int end,
    char* qname, float val) {
  // check validity of positions
  if (start < 0) {
    fprintf(stderr, "Warning! Read %s prevented from extending below 0 on %s\n",
      qname, chrom->name);
    start = 0;
  }
  if (start >= chrom->len) {
    char* msg = (char*) memalloc(MAX_ALNS);
    sprintf(msg, "Read %s, ref. %s", qname, chrom->name);
    exit(error(msg, ERRPOS));
  }
  if (end > chrom->len) {
    fprintf(stderr, "Warning! Read %s prevented from extending past %d on %s\n",
      qname, chrom->len, chrom->name);
    end = chrom->len;
  }

  // save ends of fragment to 'diff' array
  if (chrom->diff == NULL) {
    chrom->diff = (float*) memalloc((chrom->len + 1) * sizeof(float));
    for (int i = 0; i < chrom->len; i++)
      chrom->diff[i] = 0.0f;
  }
  chrom->diff[start] += val;
  chrom->diff[end] -= val;

  return end - start;
}

/* int processAvgExt()
 * Save complete intervals for unpaired alignments
 *   with "extend to average length" option, after
 *   calculating average length from paired alns.
 *   Return number printed.
 */
int processAvgExt(Aln* unpair, int unpairLen,
    double totalLen, int pairedPr) {
  // determine average fragment length
  int avgLen = 0;
  if (! pairedPr) {
    fprintf(stderr, "Warning! No paired alignments to calculate avg ");
    fprintf(stderr, "frag length --\n  Printing singletons \"as is\"\n");
  } else
    avgLen = (int) (totalLen / pairedPr + 0.5);

  int printed = 0;  // counting variable
  for (int i = 0; i < unpairLen; i++) {
    Aln* a = unpair + i;
    if (! avgLen)
      saveInterval(a->chrom, a->pos[0], a->pos[1], a->name,
        a->val);
    else if (a->strand)
      saveInterval(a->chrom, a->pos[0], a->pos[0] + avgLen,
        a->name, a->val);
    else
      saveInterval(a->chrom, a->pos[1] - avgLen, a->pos[1],
        a->name, a->val);
    printed++;

    // free memory
    free(a->name);
  }

  return printed;
}

/* void saveAvgExt()
 * Save info for an unpaired alignment to list
 *   (for "extend to average length" option), for
 *   later processing by saveAvgExt().
 */
void saveAvgExt(char* qname, Aln* b, float val,
    Aln** unpair, int* unpairLen, int* unpairMem) {

  // alloc memory if necessary
  if (*unpairLen + 1 > *unpairMem) {
    *unpairMem += MAX_ALNS;
    *unpair = (Aln*) memrealloc(*unpair,
      *unpairMem * sizeof(Aln));
  }

  // copy alignment info
  Aln* a = *unpair + *unpairLen;
  a->name = (char*) memalloc(1 + strlen(qname));
  strcpy(a->name, qname);
  a->chrom = b->chrom;
  a->strand = b->strand;
  a->pos[0] = b->pos[0];
  a->pos[1] = b->pos[1];
  a->val = val;

  (*unpairLen)++;
}

/* void saveSingle()
 * Control processing of singleton alignments
 *   (either keeping them as is, or extending
 *   to a given length).
 */
void saveSingle(char* qname, Aln* a, float val,
    bool extendOpt, int extend) {
  if (extendOpt) {
    if (a->strand)
      saveInterval(a->chrom, a->pos[0], a->pos[0] + extend,
        qname, val);
    else
      saveInterval(a->chrom, a->pos[1] - extend, a->pos[1],
        qname, val);
  } else
    saveInterval(a->chrom, a->pos[0], a->pos[1], qname, val);
}

/* int saveFragment()
 * Save full fragment for a proper pair. Return length.
 */
int saveFragment(char* qname, Aln* a, float val) {
  // ensure start < end
  int start, end;
  if (a->pos[0] > a->pos[1]) {
    start = a->pos[1];
    end = a->pos[0];
  } else {
    start = a->pos[0];
    end = a->pos[1];
  }
  return saveInterval(a->chrom, start, end, qname, val);
}

/* void processAlns()
 * Process a set of alignments. If there are properly
 *   paired alignments, process them only, weighted to
 *   1/n. Otherwise, and only if singleOpt, process
 *   singleton alignments, also weighted.
 */
void processAlns(char* qname, Aln* aln, int alnLen,
    double* totalLen, int* pairedPr, int* singlePr,
    int* orphan, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt,
    Aln** unpair, int* unpairLen, int* unpairMem) {

  // count properly paired vs. singleton alignments
  int pairCount = 0, singleR1 = 0, singleR2 = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired)
      a->full ? pairCount++ : (*orphan)++;
    else if (a->first)
      singleR1++; // singleton for R1 read
    else
      singleR2++; // singleton for R2 read
  }

  if (pairCount) {

    // save full fragments
    float val = 1.0f / pairCount; // weighted value of each aln
    int fragLen = 0;              // local sum of fragment lengths
    for (int i = 0; i < alnLen; i++) {
      Aln* a = aln + i;
      if (a->paired && a->full)
        fragLen += saveFragment(qname, a, val);
    }
    *totalLen += (double) fragLen / pairCount;
    (*pairedPr)++;

  } else if ((singleR1 || singleR2) && singleOpt) {

    // process singleton alignments (separately for R1, R2)
    for (int j = 0; j < 2; j++) {

      if ((j && ! singleR2) || (! j && ! singleR1) )
        continue;
      float val = j ? 1.0f / singleR2 : 1.0f / singleR1;  // weighted value

      // find and save singleton alignments
      for (int i = 0; i < alnLen; i++) {
        Aln* a = aln + i;
        if (! a->paired && ((j && ! a->first)
            || (! j && a->first))) {

          if (avgExtOpt) {
            // for average-extension option, save alignment
            //   for later processing by processAvgExt()
            saveAvgExt(qname, a, val, unpair,
              unpairLen, unpairMem);
          } else {
            // for other options, save singleton interval
            saveSingle(qname, a, val, extendOpt, extend);
            (*singlePr)++;
          }

        }
      }
    }
  }

}

/*** Save alignment information ***/

/* void updatePairedAln()
 * Complete a properly paired alignment.
 */
void updatePairedAln(Aln* a, uint16_t flag, uint32_t pos,
    int length) {
  if (flag & 0x40)
    a->pos[0] = flag & 0x10 ? pos + length : pos;
  else
    a->pos[1] = flag & 0x10 ? pos + length : pos;
  a->full = true;
}

/* void savePairedAln()
 * Start a properly paired alignment.
 */
void savePairedAln(Aln** aln, int* alnLen, int* alnMem,
    uint16_t flag, Chrom* chrom, uint32_t pos, int length,
    uint32_t pnext) {

  // alloc memory if necessary
  if (*alnLen + 1 > *alnMem) {
    *alnMem += MAX_ALNS;
    *aln = (Aln*) memrealloc(*aln, *alnMem * sizeof(Aln));
  }

  // save aln info
  Aln* a = *aln + *alnLen;
  a->chrom = chrom;
  a->full = false;
  a->paired = true;

  // save positions for this aln (corrected if rev-comp),
  //   and for its pair (pnext, uncorrected)
  if (flag & 0x40) {
    a->pos[0] = flag & 0x10 ? pos + length : pos;
    a->pos[1] = pnext;
  } else {
    a->pos[0] = pnext;
    a->pos[1] = flag & 0x10 ? pos + length : pos;
  }

  (*alnLen)++;
}

/* void saveSingleAln()
 * Save the information for a singleton alignment.
 */
void saveSingleAln(Aln** aln, int* alnLen, int* alnMem,
    uint16_t flag, Chrom* chrom, uint32_t pos, int length) {

  // alloc memory if necessary
  if (*alnLen + 1 > *alnMem) {
    *alnMem += MAX_ALNS;
    *aln = (Aln*) memrealloc(*aln, *alnMem * sizeof(Aln));
  }

  // save aln info
  Aln* a = *aln + *alnLen;
  a->chrom = chrom;
  a->paired = false;
  a->strand = flag & 0x10 ? false : true;
  a->first = flag & 0x40 ? true : false;
  a->pos[0] = pos;
  a->pos[1] = pos + length;
  (*alnLen)++;
}

/* void parseAlign()
 * Parse a SAM/BAM alignment record. Save alignment
 *   info to Aln* array.
 */
void parseAlign(Aln** aln, int* alnLen, int* alnMem,
    uint16_t flag, Chrom* chrom, uint32_t pos, int length,
    uint32_t pnext, int* paired, int* single, int* secPair,
    int* secSingle, bool singleOpt) {

  if ((flag & 0x3) == 0x3) {
    // paired alignment: save alignment information
    (*paired)++;
    if (flag & 0x100)
      (*secPair)++;

    // search for matching paired alignment (already analyzed)
    int i;
    for (i = 0; i < *alnLen; i++) {
      Aln* a = *aln + i;
      if (a->paired && ! a->full && a->chrom == chrom
          && (flag & 0x40 ? a->pos[0] == pos
          : a->pos[1] == pos)) {
        // complete paired alignment
        updatePairedAln(a, flag, pos, length);
        break;
      }
    }
    if (i == *alnLen)
      // start new paired alignment
      savePairedAln(aln, alnLen, alnMem, flag, chrom,
        pos, length, pnext);

  } else {

    // unpaired alignment
    (*single)++;
    if (flag & 0x100)
      (*secSingle)++;
    if (singleOpt)
      saveSingleAln(aln, alnLen, alnMem, flag, chrom,
        pos, length);

  }
}

/*** Save SAM/BAM header info ***/

/* int saveChrom()
 * If chromosome (reference sequence) has not been
 *   saved yet, save it to the array. Return the index.
 */
int saveChrom(char* name, int len, int* chromLen,
    Chrom** chrom, int xcount, char** xchrList,
    bool ctrl) {

  // determine if chrom has been saved already
  for (int i = 0; i < *chromLen; i++) {
    Chrom* c = *chrom + i;
    if (!strcmp(c->name, name)) {
      if (c->len != len)
        exit(error(c->name, ERRCHRLEN));
      return i;
    }
  }

  // determine if chrom should be skipped
  bool skip = ctrl; // automatically skip if ref in ctrl sample only
  if (! skip)
    // check if chrom is on skipped list
    for (int i = 0; i < xcount; i++)
      if (!strcmp(xchrList[i], name)) {
        skip = true;
        break;
      }

  // save to list
  *chrom = (Chrom*) memrealloc(*chrom,
    (*chromLen + 1) * sizeof(Chrom));
  Chrom* c = *chrom + *chromLen;
  c->name = (char*) memalloc(1 + strlen(name));
  strcpy(c->name, name);
  c->len = len;
  c->skip = skip;
  c->diff = NULL;
  c->treat = NULL;
  c->treatLen = 0;
  c->ctrl = NULL;
  c->ctrlLen = 0;
  c->pval = NULL;
  c->qval = NULL;
  c->pvalLen = 0;
  (*chromLen)++;
  return *chromLen - 1;
}

/* void loadChrom()
 * Save chromosome length info from a SAM header line.
 */
void loadChrom(char* line, int* chromLen, Chrom** chrom,
    int xcount, char** xchrList, bool ctrl) {
  // parse SAM header line for chrom info
  char* name = NULL, *len = NULL;
  char* field = strtok(NULL, TAB);
  while (field != NULL) {
    if (!strncmp(field, "SN:", 3))
      name = field + 3;
    else if (!strncmp(field, "LN:", 3))
      len = field + 3;
    field = strtok(NULL, TAB);
  }
  if (name == NULL || len == NULL)
    return;

  // save chrom info to array (*chrom)
  saveChrom(name, getInt(len), chromLen, chrom,
    xcount, xchrList, ctrl);
}

/* void checkHeader()
 * Check SAM header line for useful information:
 *   sort order or chromosome lengths.
 */
void checkHeader(char* line, int* chromLen, Chrom** chrom,
    int xcount, char** xchrList, bool ctrl) {

  // load tag from SAM header line
  char* tag = strtok(line, TAB);
  if (tag == NULL)
    return;

  if (! strcmp(tag, "@HD")) {
    // first header line: check sort order
    char* order = NULL;
    char* field = strtok(NULL, TAB);
    while (field != NULL) {
      if (!strncmp(field, "SO:", 3))
        order = field + 3;
      field = strtok(NULL, TAB);
    }
    // removing trailing '\n'
    int i = 0;
    for (; order[i] != '\n' && order[i] != '\0'; i++) ;
    order[i] = '\0';
    if (order == NULL || ! strcmp(order, "unknown")
        || ! strcmp(order, "coordinate"))
      exit(error("", ERRSORT));

  } else if (! strcmp(tag, "@SQ"))
    // load chrom lengths from header line
    loadChrom(line, chromLen, chrom, xcount, xchrList,
      ctrl);

}

/*** SAM parsing ***/

/* bool loadFields()
 * Load alignment info from a SAM record.
 *   Return false on failure.
 */
bool loadFields(uint16_t* flag, char** rname, uint32_t* pos,
    uint8_t* mapq, char** cigar, char** rnext, uint32_t* pnext,
    int32_t* tlen, char** seq, char** qual, char** extra) {
  int i = 2;
  char* field = strtok(NULL, TAB);
  while (field != NULL) {
    switch (i) {
      case FLAG: *flag = getInt(field); break;
      case RNAME: *rname = field; break;
      case POS: *pos = getInt(field) - 1; break;  // convert to 0-based
      case MAPQ: *mapq = getInt(field); break;
      case CIGAR: *cigar = field; break;
      case RNEXT: *rnext = field; break;
      case PNEXT: *pnext = getInt(field) - 1; break;  // convert to 0-based
      case TLEN: *tlen = getInt(field); break;
      case SEQ: *seq = field; break;
      case QUAL: *qual = field; break;
      default: return false;
      /*  *extra = (char**) memrealloc(*extra,
          (*extraCount + 1) * sizeof(char*));
        (*extra)[*extraCount] = field;
        (*extraCount)++;*/
    }
    if (++i > 11) {
      *extra = strtok(NULL, "\n");
      break;
    }
    field = strtok(NULL, TAB);
  }
  return i > 11;
}

/* int parseCigar()
 * Calculate length of sequence and offset
 *   from a CIGAR string.
 */
int parseCigar(char* cigar, int* offset) {
  int length = 0; // length of sequence
  int pos = 0;    // position of current op in cigar
  char op;
  int len = strlen(cigar);
  for (int i = 0; i < len; i++) {
    if (cigar[i] < 0x30 || cigar[i] > 0x39) {
      op = cigar[i];
      cigar[i] = '\0';
      int opLen = getInt(cigar + pos); // length of current op
      switch (op) {
        case 'M':
        case '=':
        case 'X':
          length += opLen;
          break;
        case 'I':
        case 'S':
          length += opLen;
          *offset -= opLen;
          break;
        case 'D':
          *offset += opLen;
          break;
        case 'N':
        case 'H':
        case 'P':
          break;
        default: ;
          char msg[4] = "' '";
          msg[1] = op;
          exit(error(msg, ERRCIGAR));
      }
      pos = i + 1;
    }
  }
  return length;
}

/* int calcDist()
 * Return distance to 3' end of sequence
 *   (length + offset based on CIGAR [if avl]).
 */
int calcDist(char* qname, char* seq, char* cigar) {
  int length = strcmp(seq, "*") ? strlen(seq) : 0;
  int offset = 0;
  if (strcmp(cigar, "*")) {
    int len = parseCigar(cigar, &offset);
    if (! length)
      length = len;
    else if (length != len)
      exit(error(qname, ERRMISM));
  } else if (! length)
    exit(error(qname, ERRINFO));
  return length + offset;
}

/* int readSAM()
 * Parse the alignments in a SAM file.
 */
int readSAM(File in, bool gz, char* line, Aln** aln,
    int* alnMem, char* readName, double* totalLen,
    int* unmapped, int* paired, int* single, int* pairedPr,
    int* singlePr, int* supp, int* skipped, int* lowMapQ,
    int minMapQ, int xcount, char** xchrList, int* secPair,
    int* secSingle, int* orphan, int* chromLen,
    Chrom** chrom, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, Aln** unpair,
    int* unpairMem, bool ctrl) {

  // SAM fields to save
  char* qname, *rname, *cigar, *rnext, *seq, *qual, *extra;
  uint16_t flag;
  uint32_t pos, pnext;
  int32_t tlen;
  uint8_t mapq;

  int alnLen = 0;     // number of alignments for this read
  int unpairLen = 0;  // number of singletons (with avgExtOpt)
  bool pastHeader = false;  // to check for misplaced header lines
  int count = 0;
  while (getLine(line, MAX_SIZE, in, gz) != NULL) {

    if (line[0] == '@') {
      if (pastHeader)
        exit(error(line, ERRHEAD));
      checkHeader(line, chromLen, chrom, xcount, xchrList,
        ctrl);
      continue;
    }
    pastHeader = true;

    // parse SAM record
    qname = strtok(line, TAB);
    if (qname == NULL)
      exit(error(line, ERRSAM));
    if (! loadFields(&flag, &rname, &pos, &mapq, &cigar,
        &rnext, &pnext, &tlen, &seq, &qual, &extra))
      exit(error(qname, ERRSAM));

    count++;
    if (flag & 0x4) {
      // skip unmapped
      (*unmapped)++;
      continue;
    }
    if (! strcmp(qname, "*") || ! strcmp(rname, "*")
        || pos < 0)
      // insufficient alignment info
      exit(error(qname, ERRSAM));
    if (flag & 0xE00) {
      // skip supplementary/PCR dups/low quality
      (*supp)++;
      continue;
    }
    // find matching Chrom (reference sequence)
    Chrom* ref = NULL;
    for (int i = 0; i < *chromLen; i++)
      if (! strcmp((*chrom + i)->name, rname)) {
        ref = *chrom + i;
        break;
      }
    if (ref == NULL)
      // cannot find reference sequence
      exit(error(rname, ERRCHROM));
    if (ref->skip) {
      // alignment to skipped chromosome
      (*skipped)++;
      continue;
    }
    if (mapq < minMapQ) {
      // skip low MAPQ alignments
      (*lowMapQ)++;
      continue;
    }

    // process previous set of alns, if starting a new set
    if (readName[0] == '\0' || strcmp(qname, readName)) {
      if (readName[0] != '\0')
        processAlns(readName, *aln, alnLen, totalLen,
          pairedPr, singlePr, orphan, singleOpt,
          extendOpt, extend, avgExtOpt,
          unpair, &unpairLen, unpairMem);
      alnLen = 0;
      strncpy(readName, qname, MAX_ALNS);
    }

    // save alignment information
    int length = calcDist(qname, seq, cigar); // distance to 3' end
    parseAlign(aln, &alnLen, alnMem, flag, ref, pos,
      length, pnext, paired, single, secPair, secSingle,
      singleOpt);
    // NOTE: the following SAM fields are ignored:
    //   rnext, tlen, qual, extra (optional fields)
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt,
      unpair, &unpairLen, unpairMem);

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    *singlePr += processAvgExt(*unpair, unpairLen,
      *totalLen, *pairedPr);

  return count;
}

/*** BAM parsing ***/

/* int32_t readInt32()
 * Read an int32_t in little-endian format from a
 *   gzip-compressed file. On failure, return error
 *   or EOF.
 */
int32_t readInt32(gzFile in, bool end) {
  int32_t ans = 0;
  for (int i = 0; i < sizeof(int32_t); i++) {
    int m = gzgetc(in);
    if (m == -1) {
      if (end)
        exit(error("", ERRBAM));
      else
        return EOF;
    }
    ans = ans | ((m & 0xFF) << (i*8));
  }
  return ans;
}

/* int32_t loadInt32()
 * Load an int32_t in little-endian format from a
 *   char** block. Increment *block on the fly.
 */
int32_t loadInt32(char** block) {
  int32_t ans = 0;
  for (int i = 0; i < sizeof(int32_t); i++) {
    ans = ans | ((**block & 0xFF) << (i*8));
    (*block)++;
  }
  return ans;
}

/* void loadBAMfields()
 * Load fields from a BAM record. See SAM spec section 4.2
 *   for details.
 */
void loadBAMfields(char** block, int32_t* refID, int32_t* pos,
    uint8_t* mapq, uint16_t* n_cigar_op, uint16_t* flag,
    int32_t* l_seq, int32_t* next_refID, int32_t* next_pos,
    int32_t* tlen, char** read_name, uint32_t** cigar,
    uint8_t** seq, char** qual, char** extra) {
  *refID = loadInt32(block);
  *pos = loadInt32(block);
  uint32_t bin_mq_nl = loadInt32(block);
  int8_t l_read_name = bin_mq_nl & 0xFF;
  *mapq = (bin_mq_nl >> 8) & 0xFF;
  uint32_t flag_nc = loadInt32(block);
  *n_cigar_op = flag_nc & 0xFFFF;
  *flag = (flag_nc >> 16) & 0xFFFF;
  *l_seq = loadInt32(block);
  *next_refID = loadInt32(block);
  *next_pos = loadInt32(block);
  *tlen = loadInt32(block);
  *read_name = *block;
  (*block) += l_read_name;
  *cigar = (uint32_t*) *block;
  (*block) += *n_cigar_op * sizeof(uint32_t);
  *seq = (uint8_t*) *block;
  (*block) += (*l_seq + 1)/2 * sizeof(uint8_t);
  *qual = *block;
  (*block) += *l_seq;
  *extra = *block;
}

/* int calcDistBAM()
 * Return distance to 3' end of sequence
 *   (length + offset based on BAM CIGAR).
 */
int calcDistBAM(int32_t l_seq, uint16_t n_cigar_op,
    uint32_t* cigar) {
  int length = l_seq;
  for (int i = 0; i < n_cigar_op; i++) {
    uint8_t op = cigar[i] & 0xF;
    uint32_t op_len = cigar[i] >> 4;
    if (op == 1 || op == 4)  // 'I' or 'S'
      length -= op_len;
    else if (op == 2)        // 'D'
      length += op_len;
  }
  return length;
}

/* int parseBAM()
 * Parse the alignments in a BAM file.
 */
int parseBAM(gzFile in, char* line, Aln** aln, int* alnMem,
    char* readName, int chromLen, Chrom* chrom, int n_ref,
    int idx[], double* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ, int minMapQ,
    int* secPair, int* secSingle, int* orphan,
    bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt, Aln** unpair, int* unpairMem) {

  // BAM fields to save
  int32_t refID, pos, l_seq, next_refID, next_pos, tlen;
  uint16_t n_cigar_op, flag;
  uint8_t mapq;
  uint32_t* cigar;
  uint8_t* seq;
  char* read_name, *qual, *extra;

  int alnLen = 0;     // number of alignments for this read
  int unpairLen = 0;  // number of singletons (with avgExtOpt)
  int count = 0;
  int32_t block_size;
  while ((block_size = readInt32(in, false)) != EOF) {

    if (block_size < 6 * sizeof(int32_t) + 2 * sizeof(uint32_t))
      exit(error("", ERRBAM));  // min. guaranteed block size

    // copy alignment record
    char* block = line;
    for (int i = 0; i < block_size; i++) {
      int m = gzgetc(in);
      if (m == -1)
        exit(error("", ERRBAM));
      block[i] = m;
    }

    // save BAM fields
    loadBAMfields(&block, &refID, &pos, &mapq, &n_cigar_op,
      &flag, &l_seq, &next_refID, &next_pos, &tlen,
      &read_name, &cigar, &seq, &qual, &extra);
    if (block > line + block_size)
      exit(error("", ERRBAM));  // read off the end of the block

    count++;
    if (flag & 0x4) {
      // skip unmapped
      (*unmapped)++;
      continue;
    }
    if (! strcmp(read_name, "*")
        || refID < 0 || refID >= n_ref
        || idx[refID] < 0 || idx[refID] >= chromLen
        || pos < 0)
      // insufficient alignment info
      exit(error(read_name, ERRSAM));
    if (flag & 0xE00) {
      // skip supplementary/PCR dups/low quality
      (*supp)++;
      continue;
    }
    Chrom* ref = chrom + idx[refID];
    if (ref->skip) {
      // alignment to skipped chromosome
      (*skipped)++;
      continue;
    }
    if (mapq < minMapQ) {
      // skip low MAPQ alignments
      (*lowMapQ)++;
      continue;
    }

    // process previous set of alns
    if (readName[0] == '\0' || strcmp(read_name, readName)) {
      if (readName[0] != '\0')
        processAlns(readName, *aln, alnLen, totalLen,
          pairedPr, singlePr, orphan, singleOpt,
          extendOpt, extend, avgExtOpt,
          unpair, &unpairLen, unpairMem);
      alnLen = 0;
      strncpy(readName, read_name, MAX_ALNS);
    }

    // save alignment information
    int length = calcDistBAM(l_seq, n_cigar_op, cigar); // distance to 3' end
    parseAlign(aln, &alnLen, alnMem, flag, ref, pos,
      length, next_pos, paired, single, secPair, secSingle,
      singleOpt);
    // NOTE: the following BAM fields are ignored:
    //   next_refID, tlen, seq, qual, extra (optional fields)
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt,
      unpair, &unpairLen, unpairMem);

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    *singlePr += processAvgExt(*unpair, unpairLen,
      *totalLen, *pairedPr);

  return count;
}

/* int readBAM()
 * Parse the header from a BAM file, then
 *   call parseBAM().
 */
int readBAM(gzFile in, char* line, Aln** aln, int* alnMem,
    char* readName, double* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ, int minMapQ,
    int xcount, char** xchrList, int* secPair,
    int* secSingle, int* orphan, int* chromLen,
    Chrom** chrom, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, Aln** unpair,
    int* unpairMem, bool ctrl) {

  // load first line from header
  int32_t l_text = readInt32(in, true);
  int i;
  for (i = 0; i < l_text; i++) {
    int m = gzgetc(in);
    if (m == -1)
      exit(error("", ERRBAM));
    unsigned char n = m;
    if (n == '\n' || n == '\0')
      break;
    line[i] = n;
  }
  line[i] = '\0';
  // check sort order
  char* tag = strtok(line, TAB);
  if (tag == NULL || strcmp(tag, "@HD"))
    exit(error("", ERRBAM));
  char* order = NULL;
  char* field = strtok(NULL, TAB);
  while (field != NULL) {
    if (!strncmp(field, "SO:", 3))
      order = field + 3;
    field = strtok(NULL, TAB);
  }
  if (order == NULL || ! strcmp(order, "unknown")
      || ! strcmp(order, "coordinate"))
    exit(error("", ERRSORT));
  if (gzseek(in, l_text - i - 1, SEEK_CUR) == -1)
    exit(error("", ERRBAM));

  // save chromosome lengths
  int32_t n_ref = readInt32(in, true); // number of ref sequences
  int idx[n_ref];  // index of reference sequences into *chrom array
  for (int i = 0; i < n_ref; i++) {
    int32_t len = readInt32(in, true);
    if (len < 1 || len > MAX_SIZE)
      exit(error("", ERRBAM));
    for (int j = 0; j < len; j++) {
      int m = gzgetc(in);
      if (m == -1)
        exit(error("", ERRBAM));
      line[j] = m;
    }
    if (line[len-1] != '\0')
      exit(error("", ERRBAM));
    idx[i] = saveChrom(line, readInt32(in, true),
      chromLen, chrom, xcount, xchrList, ctrl);
  }

  return parseBAM(in, line, aln, alnMem, readName,
    *chromLen, *chrom, n_ref, idx, totalLen, unmapped,
    paired, single, pairedPr, singlePr, supp, skipped,
    lowMapQ, minMapQ, secPair, secSingle, orphan,
    singleOpt, extendOpt, extend, avgExtOpt,
    unpair, unpairMem);
}

/*** File I/O ***/

/* void openWrite()
 * Open a file for writing (stdout if file is '-').
 */
void openWrite(char* outFile, File* out, bool gz) {
  if (outFile[0] == '-' && strlen(outFile) > 1)
    exit(error(outFile, ERRNAME));
  if (gz) {
    if (!strcmp(outFile + strlen(outFile) - strlen(GZEXT), GZEXT)
        || !strcmp(outFile, "/dev/null"))
      out->gzf = gzopen(outFile, "w");
    else if (!strcmp(outFile, "-"))
      out->gzf = gzdopen(fileno(stdout), "wb");
    else {
      // add ".gz" to outFile
      char* outFile2 = memalloc(strlen(outFile)
        + strlen(GZEXT) + 1);
      strcpy(outFile2, outFile);
      strcat(outFile2, GZEXT);
      out->gzf = gzopen(outFile2, "w");
      free(outFile2);
    }
    if (out->gzf == NULL)
      exit(error(outFile, ERROPENW));
  } else {
    out->f = (strcmp(outFile, "-") ?
      fopen(outFile, "w") : stdout);
    if (out->f == NULL)
      exit(error(outFile, ERROPENW));
  }
}

/* void openFiles()
 * Opens output files for the program,
 *   adjusting file names/extensions as needed.
 */
void openFiles(char* outFile, File* out,
    char* logFile, File* log, bool gz, bool qvalOpt) {

  // open peak file
  openWrite(outFile, out, gz);

  // open bedgraph file
  if (logFile != NULL) {
    openWrite(logFile, log, gz);

    // print header
    if (gz) {
      gzprintf(log->gzf, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)");
      if (qvalOpt)
        gzprintf(log->gzf, "\t-log(q)");
      gzprintf(log->gzf, "\tsignif\n");
    } else {
      fprintf(log->f, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)");
      if (qvalOpt)
        fprintf(log->f, "\t-log(q)");
      fprintf(log->f, "\tsignif\n");
    }
  }

}

/* bool checkBAM()
 * Determine if file is BAM formatted.
 */
bool checkBAM(File in, bool gz) {
  if (! gz)
    return false;
  char magic[5] = "BAM\1";  // BAM magic string
  for (int i = 0; i < 4; i++) {
    char m = gzgetc(in.gzf);
    if (m == -1)
      exit(error("", ERROPEN));
    if (m != magic[i]) {
      // push back chars before returning false
      if (gzungetc(m, in.gzf) == -1)
        exit(error("", ERRUNGET));
      for (int j = i - 1; j > -1; j--)
        if (gzungetc(magic[j], in.gzf) == -1)
          exit(error("", ERRUNGET));
      return false;
    }
  }
  return true;
}

/* bool openRead()
 * Open a file for reading (stdin if file is '-').
 *   Return true if gzip compressed.
 */
bool openRead(char* inFile, File* in) {

  // open file or stdin
  bool stdinBool = (strcmp(inFile, "-") ? false : true);
  FILE* dummy = (stdinBool ? stdin : fopen(inFile, "r"));
  if (dummy == NULL)
    exit(error(inFile, ERROPEN));

  // check for gzip compression: magic number 0x1F, 0x8B
  bool gzip = true;
  int save = 0;  // first char to pushback (for stdin)
  int i, j;
  for (i = 0; i < 2; i++) {
    j = fgetc(dummy);
    if (j == EOF)
      exit(error(inFile, ERROPEN));
    if ( (i && (unsigned char) j != 0x8B)
        || (! i && (unsigned char) j != 0x1F) ) {
      gzip = false;
      break;
    }
    if (! i)
      save = j;
  }

  // for stdin, push back chars
  if (stdinBool) {
    if (gzip)
      exit(error("", ERRGZIP));
    if (ungetc(j, dummy) == EOF)
      exit(error("", ERRUNGET));
    if (i && ungetc(save, dummy) == EOF)
      exit(error("", ERRUNGET));
  }

  // open file
  if (gzip) {
    if (fclose(dummy))
      exit(error("<dummy>", ERRCLOSE));
    in->gzf = gzopen(inFile, "r");
    if (in->gzf == NULL)
      exit(error(inFile, ERROPEN));
  } else {
    if (! stdinBool)
      rewind(dummy);
    in->f = dummy;
  }

  return gzip;
}

/* void logCounts()
 * Log alignment counts to stderr.
 */
void logCounts(int count, int unmapped, int supp,
    int skipped, Chrom* chrom, int chromLen, int minMapQ,
    int lowMapQ, int paired, int secPair, int orphan,
    int single, int secSingle, int singlePr, int pairedPr,
    double totalLen, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, bool bam) {
  double avgLen = pairedPr ? totalLen / pairedPr : 0.0;
  fprintf(stderr, "  %s records analyzed: %10d\n",
    bam ? "BAM" : "SAM", count);
  if (unmapped)
    fprintf(stderr, "    Unmapped:           %10d\n", unmapped);
  if (supp)
    fprintf(stderr, "    Supp./dups/lowQual: %10d\n", supp);
  if (skipped) {
    fprintf(stderr, "    To skipped refs:    %10d\n", skipped);
    fprintf(stderr, "      (");
    bool first = true;
    for (int i = 0; i < chromLen; i++) {
      Chrom* c = chrom + i;
      if (c->skip) {
        fprintf(stderr, "%s%s", first ? "" : ",", c->name);
        first = false;
      }
    }
    fprintf(stderr, ")\n");
  }
  if (lowMapQ)
    fprintf(stderr, "    MAPQ < %-2d:          %10d\n", minMapQ, lowMapQ);
  fprintf(stderr, "    Paired alignments:  %10d\n", paired);
  if (secPair)
    fprintf(stderr, "      secondary alns:   %10d\n", secPair);
  if (orphan)
    fprintf(stderr, "      \"orphan\" alns:    %10d\n", orphan);
  fprintf(stderr, "    Unpaired alignments:%10d\n", single);
  if (secSingle)
    fprintf(stderr, "      secondary alns:   %10d\n", secSingle);
  fprintf(stderr, "  Fragments analyzed:   %10d\n", singlePr + pairedPr);
  fprintf(stderr, "    Full fragments:     %10d\n", pairedPr);
  if (pairedPr)
    fprintf(stderr, "      (avg. length: %.1fbp)\n", avgLen);
  if (singleOpt) {
    fprintf(stderr, "    Singletons:         %10d\n", singlePr);
    if (singlePr) {
      if (extendOpt)
        fprintf(stderr, "      (extended to length %dbp)\n", extend);
      else if (avgExtOpt && pairedPr)
        fprintf(stderr, "      (extended to length %dbp)\n",
          (int) (avgLen + 0.5));
    }
  }
}

/* void runProgram()
 * Controls the opening/closing of files,
 *   and analysis by readSAM() or readBAM().
 *   Passes results to findPeaks().
 */
void runProgram(char* inFile, char* ctrlFile, char* outFile,
    char* logFile, bool gzOut, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, int minMapQ, int xcount,
    char** xchrList, float pqvalue, bool qvalOpt,
    int minLen, int maxGap, bool verbose, int threads) {

  // initialize variables
  char* line = (char*) memalloc(MAX_SIZE);
  int chromLen = 0;       // number of reference sequences
  Chrom* chrom = NULL;    // array of reference sequences
  int alnMem = MAX_ALNS;  // number of alignments to save per read (dynamic)
  Aln* aln = (Aln*) memalloc(alnMem * sizeof(Aln)); // array of saved alns
  int unpairMem = 0;      // number of unpaired alns (for avg-ext option)
  Aln* unpair = NULL;     // array of unpaired alns (for avg-ext option)
  char* readName = memalloc(MAX_ALNS + 1);  // name of read being analyzed
  readName[0] = readName[MAX_ALNS] = '\0';
  double fragLen = 0.0;   // total weighted length of all treatment fragments
  unsigned long genomeLen = 0;  // total length of genome

  // loop through input files (treatment and control)
  for (int i = 0; i < 2; i++) {

    // get treat/ctrl file name(s)
    char* file = (i ? ctrlFile : inFile);
    if (file == NULL) {
      if (verbose)
        fprintf(stderr, "No control file provided\n");
      savePileupNoCtrl(chrom, chromLen, fragLen, &genomeLen);
      break;
    }

    char* end;
    char* filename = strtok_r(file, COM, &end);
    while (filename) {
      // open input file
      File in;
      bool gz = openRead(filename, &in);
      bool bam = checkBAM(in, gz);

      // process files
      if (verbose)
        fprintf(stderr, "Processing %s file: %s\n",
          i ? "control" : "treatment", filename);
      int unmapped = 0, paired = 0, single = 0, orphan = 0,
        pairedPr = 0, singlePr = 0, supp = 0, skipped = 0,
        lowMapQ = 0, secPair = 0, secSingle = 0;  // counting variables
      double totalLen = 0.0; // total weighted length of paired fragments
      int count;
      if (bam)
        count = readBAM(in.gzf, line, &aln, &alnMem,
          readName, &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ,
          minMapQ, xcount, xchrList, &secPair, &secSingle,
          &orphan, &chromLen, &chrom, singleOpt, extendOpt,
          extend, avgExtOpt, &unpair, &unpairMem, i);
      else
        count = readSAM(in, gz, line, &aln, &alnMem,
          readName, &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ,
          minMapQ, xcount, xchrList, &secPair, &secSingle,
          &orphan, &chromLen, &chrom, singleOpt, extendOpt,
          extend, avgExtOpt, &unpair, &unpairMem, i);

      // log counts
      if (verbose)
        logCounts(count, unmapped, supp, skipped, chrom,
          chromLen, minMapQ, lowMapQ, paired, secPair,
          orphan, single, secSingle, singlePr, pairedPr,
          totalLen, singleOpt, extendOpt, extend,
          avgExtOpt, bam);

      // close input files
      if ( (gz && gzclose(in.gzf) != Z_OK) || (! gz && fclose(in.f)) )
        exit(error(filename, ERRCLOSE));

      filename = strtok_r(NULL, COM, &end);
    }

    // save pileup values
    if (i)
      savePileupCtrl(chrom, chromLen, fragLen, &genomeLen,
        1.0f / alnMem);
    else {
      fragLen = savePileupTreat(chrom, chromLen,
        1.0f / alnMem);
      // reset 'diff' array for each Chrom
      for (int j = 0; j < chromLen; j++) {
        Chrom* c = chrom + j;
        if (c->diff != NULL)
          for (int k = 0; k < c->len + 1; k++)
            c->diff[k] = 0.0f;
      }
    }
  }

  // open output files
  File out, log;
  openFiles(outFile, &out, logFile, &log, gzOut, qvalOpt);

  // find peaks
  findPeaks(out, log, logFile != NULL, gzOut, genomeLen,
    chrom, chromLen, pqvalue, qvalOpt, minLen, maxGap,
    verbose);

  // free memory
  if (xcount) {
    for (int i = 0; i < xcount; i++)
      free(xchrList[i]);
    free(xchrList);
  }
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip) {
      if (qvalOpt) {
        free(chr->qval->end);
        free(chr->qval->cov);
        free(chr->qval);
      }
      free(chr->pval->end);
      free(chr->pval->cov);
      free(chr->pval);
      free(chr->ctrl->end);
      free(chr->ctrl->cov);
      free(chr->ctrl);
      free(chr->treat->end);
      free(chr->treat->cov);
      free(chr->treat);
      free(chr->diff);
    }
    free(chr->name);
  }
  free(chrom);
  free(aln);
  free(unpair);
  free(readName);
  free(line);

  // close files
  if ( ( gzOut && gzclose(out.gzf) != Z_OK ) ||
      ( ! gzOut && fclose(out.f) ) )
    exit(error(outFile, ERRCLOSE));
  if (logFile != NULL && ( ( gzOut && gzclose(log.gzf) != Z_OK )
      || ( ! gzOut && fclose(log.f) ) ) )
    exit(error(logFile, ERRCLOSE));
}

/* int saveXChrom()
 * Save list of chromosomes (ref names) to ignore.
 *   Return count.
 */
int saveXChrom(char* xchrom, char*** xchrList) {
  int i = 0;
  char* chrom = strtok(xchrom, COM);
  while (chrom != NULL) {
    *xchrList = (char**) memrealloc(*xchrList, (i+1) * sizeof(char*));
    (*xchrList)[i] = (char*) memalloc(1 + strlen(chrom));
    strcpy((*xchrList)[i], chrom);
    i++;
    chrom = strtok(NULL, COM);
  }
  return i;
}

/* void getArgs()
 * Parse the command-line. Check for errors.
 */
void getArgs(int argc, char** argv) {

  // default parameters/filenames
  char* outFile = NULL, *inFile = NULL,
    *ctrlFile = NULL, *logFile = NULL;
  char* xchrom = NULL;
  int extend = 0, minMapQ = 0, minLen = DEFMINLEN,
    maxGap = DEFMAXGAP, threads = DEFTHR;
  float pqvalue = DEFQVAL;
  bool singleOpt = false, extendOpt = false,
    avgExtOpt = false, gzOut = false, qvalOpt = true;
  bool verbose = false;

  // parse argv
  int c;
  while ( (c = getopt_long(argc, argv, OPTIONS, long_options, NULL)) != -1 )
    switch (c) {
      case INFILE: inFile = optarg; break;
      case CTRLFILE: ctrlFile = optarg; break;
      case OUTFILE: outFile = optarg; break;
      case LOGFILE: logFile = optarg; break;
      case GZOPT: gzOut = true; break;
      case SINGLEOPT: singleOpt = true; break;
      case EXTENDOPT: extend = getInt(optarg); extendOpt = true; break;
      case AVGEXTOPT: avgExtOpt = true; break;
      case XCHROM: xchrom = optarg; break;
      case MINMAPQ: minMapQ = getInt(optarg); break;
      case QVALUE: pqvalue = getFloat(optarg); break;
      case PVALUE: pqvalue = getFloat(optarg); qvalOpt = false; break;
      case MINLEN: minLen = getInt(optarg); break;
      case MAXGAP: maxGap = getInt(optarg); break;

      case VERBOSE: verbose = true; break;
      case VERSOPT: printVersion(); break;
      case HELP: usage(); break;

      case THREADS: threads = getInt(optarg); break;
      default: exit(-1);
    }
  if (optind < argc)
    exit(error(argv[optind], ERRPARAM));

  // check for argument errors
  if (outFile == NULL || inFile == NULL) {
    error("", ERRFILE);
    usage();
  }
  if (avgExtOpt) {
    singleOpt = true;
    extendOpt = false; // avgExtOpt takes precedence
  }
  if (extendOpt) {
    singleOpt = true;
    if (extend <= 0)
      exit(error("", ERREXTEND));
  }
  if (threads < 1)
    exit(error("", ERRTHREAD));

  // save list of chromosomes to ignore
  int xcount = 0;
  char** xchrList = NULL;
  if (xchrom != NULL)
    xcount = saveXChrom(xchrom, &xchrList);

  // adjust significance level to -log scale
  pqvalue = -1 * log10f(pqvalue);

  // send arguments to runProgram()
  runProgram(inFile, ctrlFile, outFile, logFile, gzOut,
    singleOpt, extendOpt, extend, avgExtOpt, minMapQ,
    xcount, xchrList, pqvalue, qvalOpt, minLen, maxGap,
    verbose, threads);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getArgs(argc, argv);
  return 0;
}
