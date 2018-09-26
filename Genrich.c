/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.3
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <limits.h>
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
  fprintf(stderr, "                     (matched with -%c files; 'null' if missing)\n", INFILE);
  fprintf(stderr, "  -%c  <file>       Output bedgraph-ish file for p/q values\n", LOGFILE);
  fprintf(stderr, "  -%c  <file>       Output bedgraph-ish file for pileups and p-values\n", PILEFILE);
  fprintf(stderr, "  -%c  <file>       Output BED file for reads/fragments/intervals\n", BEDFILE);
  fprintf(stderr, "Filtering options:\n");
  fprintf(stderr, "  -%c  <arg>        Comma-separated list of chromosomes to ignore\n", XCHROM);
  fprintf(stderr, "  -%c  <file>       Input BED file of genomic regions to ignore\n", XFILE);
  fprintf(stderr, "  -%c  <int>        Minimum MAPQ to keep an alignment (def. 0)\n", MINMAPQ);
  fprintf(stderr, "  -%c  <float>      Keep secondary alignments whose scores are\n", ASDIFF);
  fprintf(stderr, "                     within <float> of the best (def. 0)\n");
  fprintf(stderr, "Options for unpaired alignments:\n");
  fprintf(stderr, "  -%c               Keep unpaired alignments (def. false)\n", SINGLEOPT);
  fprintf(stderr, "  -%c  <int>        Keep unpaired alignments, with fragment length\n", EXTENDOPT);
  fprintf(stderr, "                     increased to specified value\n");
  fprintf(stderr, "  -%c               Keep unpaired alignments, with fragment length\n", AVGEXTOPT);
  fprintf(stderr, "                     increased to average value of paired alignments\n");
  fprintf(stderr, "Options for ATAC-seq:\n");
  fprintf(stderr, "  -%c               Use ATAC-seq mode (def. false)\n", ATACOPT);
  fprintf(stderr, "  -%c  <int>        Create intervals of this length, in bp (def. %d)\n", ATACLEN, DEFATAC);
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
int error(const char* msg, enum errCode err) {
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
  char* endptr;
  float ans = strtof(in, &endptr);
  if (*endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char* endptr;
  int ans = (int) strtol(in, &endptr, 10);
  if (*endptr != '\0')
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

/*** Quicksort (of p-values, for q-value calculation) ***/
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
void swapInt(uint64_t* a, uint64_t* b) {
  uint64_t t = *a;
  *a = *b;
  *b = t;
}
int64_t partition(float* pVal, uint64_t* pEnd,
    int64_t low, int64_t high) {
  float pivot = pVal[high];  // pivot value: last elt
  int64_t idx = low - 1;

  for (int64_t j = low; j < high; j++) {
    if (pVal[j] < pivot) {
      idx++;
      swapFloat(pVal + idx, pVal + j);
      swapInt(pEnd + idx, pEnd + j);  // swap int values too
    }
  }
  idx++;
  swapFloat(pVal + idx, pVal + high);
  swapInt(pEnd + idx, pEnd + high);
  return idx;
}
void quickSort(float* pVal, uint64_t* pEnd,
    int64_t low, int64_t high) {
  if (low < high) {
    int64_t idx = partition(pVal, pEnd, low, high);
    quickSort(pVal, pEnd, low, idx - 1);
    quickSort(pVal, pEnd, idx + 1, high);
  }
}

/*** Calculate q-values ***/

/* float lookup()
 * Return the pre-computed q-value for a given p-value,
 *   using parallel arrays (pVal and qVal).
 */
float lookup(float* pVal, uint64_t low, uint64_t high,
    float* qVal, float p) {
  if (low == high)
    return qVal[low];
  uint64_t idx = (low + high) / 2;
  if (pVal[idx] == p)
    return qVal[idx];
  if (pVal[idx] > p)
    return lookup(pVal, low, idx - 1, qVal, p);
  return lookup(pVal, idx + 1, high, qVal, p);
}

/* void saveQval()
 * Calculate and save q-values, given the pre-compiled
 *   arrays of p-values (pVal) and lengths (pEnd).
 */
void saveQval(Chrom* chrom, int chromLen, int n,
    uint64_t genomeLen, float* pVal,
    uint64_t* pEnd, int64_t pLen) {

  // sort pileup by p-values
  quickSort(pVal, pEnd, 0, pLen - 1);

  // calculate q-values for each p-value: -log(q) = -log(p*N/k)
  uint64_t k = 1;  // 1 + number of bases with higher -log(p)
  float logN = - log10f(genomeLen);
  float* qVal = (float*) memalloc((pLen + 1) * sizeof(float));
  qVal[pLen] = FLT_MAX;
  for (int64_t i = pLen - 1; i > -1; i--) {
    // ensure monotonicity
    qVal[i] = MAX( MIN( pVal[i] + logN + log10f(k),
      qVal[i + 1]), 0.0f);
    k += pEnd[i];
  }

  // save pileups of q-values for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || chr->pval[n] == NULL)
      continue;
    for (uint32_t j = 0; j < chr->pvalLen[n]; j++)
      chr->qval->cov[j] = lookup(pVal, 0, pLen, qVal,
        chr->pval[n]->cov[j]);
  }

  // free memory
  free(qVal);
}

/*** Save p-values in hashtable ***/

/* uint32_t jenkins_one_at_a_time_hash()
 * Adapted from http://www.burtleburtle.net/bob/hash/doobs.html
 *   Modified to take a float (p-value) as input.
 *   Returns index into hashtable.
 */
uint32_t jenkins_one_at_a_time_hash(float f) {
  uint32_t hash = 0;
  unsigned char* p = (unsigned char*) &f;
  for (int i = 0; i < sizeof(float); i++) {
    hash += p[i];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;
  return hash % HASH_SIZE;
}

/* int recordPval()
 * Save length of given p-value into hashtable.
 *   Return 1 if new entry made, else 0.
 */
int recordPval(Hash** table, float p, uint32_t length) {

  // check hashtable for matching p-value
  uint32_t idx = jenkins_one_at_a_time_hash(p);
  for (Hash* h = table[idx]; h != NULL; h = h->next)
    if (p == h->val) {
      // match: add length to bucket
      h->len += length;
      return 0;
    }

  // no match: add info into bucket
  Hash* newVal = (Hash*) memalloc(sizeof(Hash));
  newVal->val = p;
  newVal->len = length;
  newVal->next = table[idx];
  table[idx] = newVal;
  return 1;
}

/* float* collectPval()
 * Collect arrays of p-values and genome lengths from
 *   hashtable (to be used in q-value calculations).
 */
float* collectPval(Hash** table, uint64_t** pEnd,
    int64_t pLen) {
  float* pVal = (float*) memalloc(pLen * sizeof(float));
  int64_t idx = 0;
  for (int i = 0; i < HASH_SIZE; i++)
    for (Hash* h = table[i]; h != NULL; h = h->next) {
      pVal[idx] = h->val;
      (*pEnd)[idx] = h->len;
      idx++;
    }
  if (idx != pLen)
    exit(error(errMsg[ERRPVAL], ERRISSUE));
  return pVal;
}

/* Hash** hashPval()
 * Collect p-values in a hashtable.
 */
Hash** hashPval(Chrom* chrom, int chromLen, int n,
    int64_t* pLen) {

  // create hashtable for conversion of p-values to q-values
  Hash** table = (Hash**) memalloc(HASH_SIZE * sizeof(Hash*));
  for (int i = 0; i < HASH_SIZE; i++)
    table[i] = NULL;

  // loop through chroms
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || chr->pval[n] == NULL)
      continue;

    // populate hashtable
    Pileup* p = chr->pval[n]; // use the last p-value array
    uint32_t start = 0;
    for (uint32_t m = 0; m < chr->pvalLen[n]; m++) {
      // record p-value and length in hashtable
      *pLen += recordPval(table, p->cov[m],
        p->end[m] - start);
      start = p->end[m];
    }

  }

  return table;
}

/* void computeQval()
 * Control q-value calculations.
 */
void computeQval(Chrom* chrom, int chromLen,
    uint64_t genomeLen, int n) {

  // create "pileup" arrays for q-values
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || chr->pval[n] == NULL)
      continue;
    uint32_t num = chr->pvalLen[n]; // use last p-value array length
    chr->qval = (Pileup*) memalloc(sizeof(Pileup));
    chr->qval->end = (uint32_t*) memalloc(num * sizeof(uint32_t));
    chr->qval->cov = (float*) memalloc(num * sizeof(float));
  }

  // save all p-values (genome-wide) to hashtable
  int64_t pLen = 0;
  Hash** table = hashPval(chrom, chromLen, n, &pLen);

  // collect p-values from hashtable
  uint64_t* pEnd = memalloc(pLen * sizeof(uint64_t));
  float* pVal = collectPval(table, &pEnd, pLen);

  // convert p-values to q-values
  saveQval(chrom, chromLen, n, genomeLen, pVal, pEnd, pLen);

  // free memory
  free(pEnd);
  free(pVal);
  for (int i = 0; i < HASH_SIZE; i++) {
    Hash* tmp;
    Hash* h = table[i];
    while (h != NULL) {
      tmp = h->next;
      free(h);
      h = tmp;
    }
  }
  free(table);
}

/*** Calculate p-value for Chi-squared test ***/
// adapted from R-3.5.0 source code, as noted below

// from dpq.h in R-3.5.0:
#define R_Log1_Exp(x)  ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))

/* double bd0()
 * Adapted from bd0.c in R-3.5.0.
 */
double bd0(double x, double np) {
  double ej, s, s1, v;
  if (fabs(x-np) < 0.1*(x+np)) {
    v = (x-np)/(x+np);
    s = (x-np)*v;
    if (fabs(s) < DBL_MIN)
      return s;
    ej = 2*x*v;
    v = v*v;
    for (int j = 1; j < 1000; j++) {
      ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1 == s)
        return s1;
      s = s1;
    }
  }
  return x * log(x / np) + np - x;
}

/* double stirlerr()
 * Adapted from stirlerr.c in R-3.5.0.
 *   Argument 'n' is an integer, 1 <= n <= 199.
 */
double stirlerr(double n) {
  double S0 = 1.0 / 12;
  double S1 = 1.0 / 360;
  double S2 = 1.0 / 1260;
  double S3 = 1.0 / 1680;
  double S4 = 1.0 / 1188;
  double sferr[16] = {
    0.0,
    0.0810614667953272582196702,
    0.0413406959554092940938221,
    0.02767792568499833914878929,
    0.02079067210376509311152277,
    0.01664469118982119216319487,
    0.01387612882307074799874573,
    0.01189670994589177009505572,
    0.010411265261972096497478567,
    0.009255462182712732917728637,
    0.008330563433362871256469318,
    0.007573675487951840794972024,
    0.006942840107209529865664152,
    0.006408994188004207068439631,
    0.005951370112758847735624416,
    0.005554733551962801371038690
  };

  double nn = n * n;
  if (n > 80.0)
    return (S0-(S1-S2/nn)/nn)/n;
  if (n > 35.0)
    return (S0-(S1-(S2-S3/nn)/nn)/nn)/n;
  if (n > 15.0)
    return (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n;
  return sferr[(int) n];
}

/* double dpois()
 * Adapted from dpois.c in R-3.5.0 (cf. dpois_raw()).
 */
double dpois(double x, double lambda) {
  return -0.5 * log(2.0 * M_PI * x) - stirlerr(x)
    - bd0(x, lambda);
}

/* double pd_upper_series()
 * Adapted from pgamma.c in R-3.5.0.
 */
double pd_upper_series(double x, double alph) {
  double term = x / alph;
  double sum = term;
  do {
    alph++;
    term *= x / alph;
    sum += term;
  } while (term > sum * DBL_EPSILON);
  return log(sum);
}

/* double pd_lower_series()
 * Adapted from pgamma.c in R-3.5.0.
 */
double pd_lower_series(double lambda, double y) {
  double term = 1, sum = 0;
  while (y >= 1 && term > sum * DBL_EPSILON) {
    term *= y / lambda;
    sum += term;
    y--;
  }
  return log1p(sum);
}

/* double pgamma_smallx()
 * Adapted from pgamma.c in R-3.5.0.
 */
double pgamma_smallx(double x, double alph) {
  double sum = 0.0;
  double c = alph;
  double n = 0.0;
  double term;
  do {
    n++;
    c *= -x / n;
    term = c / (alph + n);
    sum += term;
  } while (fabs(term) > DBL_EPSILON * fabs(sum));
  double lf2 = alph * log(x) - lgamma(alph + 1);
  return R_Log1_Exp(log1p(sum) + lf2);
}

/* double pgamma()
 * Adapted from pgamma.c in R-3.5.0 (cf. pgamma_raw()).
 *   Argument 'alph' is an integer, 2 <= alph <= 200.
 */
double pgamma(double x, double alph) {

  if (x < 1)
    // small values of x
    return pgamma_smallx(x, alph);

  else if (x <= alph - 1) {
    // larger alph than x
    double sum = pd_upper_series(x, alph);
    double d = dpois(alph - 1, x);
    return R_Log1_Exp(sum + d);
  }

  // x > alph - 1
  double sum = pd_lower_series(x, alph - 1);
  double d = dpois(alph - 1, x);
  return sum + d;
}

/* double pchisq()
 * Calculate a p-value for a chi-squared distribution
 *   with observation 'x' and 'df' degrees of freedom.
 *   'df' must be an even integer, 4 <= df <= 400.
 * Adapted from pchisq.c and pgamma.c in R-3.5.0,
 *   with lower_tail=FALSE and log_p=TRUE.
 * Return value is -log10(p).
 */
double pchisq(double x, int df) {
  if (df < 4 || df > 400 || df / 2.0 != (int) (df / 2.0))
    exit(error(errMsg[ERRDF], ERRISSUE));
  return -pgamma(x / 2.0, df / 2.0) / M_LN10;
}

/*** Combine p-values from multiple replicates ***/

/* float multPval()
 * Combine multiple p-values into a single net p-value
 *   using Fisher's method.
 */
float multPval(Pileup** pval, int n, uint32_t idx[]) {
  double sum = 0.0;
  int df = 0;
  for (int j = 0; j < n; j++)
    if (pval[j] != NULL) {
      sum += pval[j]->cov[idx[j]];
      df += 2;
    }
  if (! df)
    exit(error(errMsg[ERRMULT], ERRISSUE));
  if (df == 2 || ! sum)
    return (float) sum;

  // calculate p-value using chi-squared dist.
  double p = pchisq(2.0 * sum / M_LOG10E, df);
  return p > FLT_MAX ? FLT_MAX : (float) p;
}

/* uint32_t countIntervals2()
 * Count the number of pileup intervals to create
 *   for the combined p-values.
 */
uint32_t countIntervals2(Chrom* c, int n) {
  uint32_t num = 1;
  uint32_t idx[n];  // indexes into each pval array
  for (int j = 0; j < n; j++)
    idx[j] = 0;
  for (uint32_t k = 1; k < c->len; k++) {
    bool add = false;
    for (int j = 0; j < n; j++)
      if (c->pval[j] != NULL
          && c->pval[j]->end[idx[j]] == k) {
        if (! add) {
          num++;
          add = true;
        }
        idx[j]++;
      }
  }
  return num;
}

/* void combinePval()
 * Combine p-values for multiple replicates.
 */
void combinePval(Chrom* chrom, int chromLen, int n) {

  // combine p-value "pileups" for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // make sure at least one pval array exists
    int j;
    for (j = 0; j < n; j++)
      if (chr->pval[j] != NULL)
        break;
    if (j == n) {
      // none exists: append another NULL
      chr->pval = (Pileup**) memrealloc(chr->pval,
        (n + 1) * sizeof(Pileup*));
      chr->pval[n] = NULL;
      continue;
    }

    // create additional 'pileup' array for combined p-values
    uint32_t num = countIntervals2(chr, n);
    chr->pval = (Pileup**) memrealloc(chr->pval,
      (n + 1) * sizeof(Pileup*));
    chr->pval[n] = (Pileup*) memalloc(sizeof(Pileup));
    chr->pval[n]->end = (uint32_t*) memalloc(num * sizeof(uint32_t));
    chr->pval[n]->cov = (float*) memalloc(num * sizeof(float));
    chr->pvalLen = (uint32_t*) memrealloc(chr->pvalLen,
      (n + 1) * sizeof(uint32_t));
    chr->pvalLen[n] = num;
    chr->sample++;

    // save combined p-values
    uint32_t idx[n + 1];  // indexes into each pval array
    for (int j = 0; j <= n; j++)
      idx[j] = 0;
    for (uint32_t k = 1; k <= chr->len; k++) {
      bool add = false;
      for (int j = 0; j < n; j++)
        if (chr->pval[j] != NULL
            && chr->pval[j]->end[idx[j]] == k) {
          if (! add) {
            chr->pval[n]->end[idx[n]] = k;
            chr->pval[n]->cov[idx[n]]
              = multPval(chr->pval, n, idx);
            idx[n]++;
            add = true;
          }
          idx[j]++;
        }

    }

  }
}

/*** Call peaks ***/

/* void printLogHeader()
 * Print header of logfile.
 */
void printLogHeader(File log, bool gzOut, int n,
    bool qvalOpt) {
  if (n) {
    // multiple samples: logfile has multiple p-values, no pileups
    if (gzOut) {
      gzprintf(log.gzf, "chr\tstart\tend");
      for (int i = 0; i < n; i++)
        gzprintf(log.gzf, "\t-log(p)_%d", i);
      gzprintf(log.gzf, "\t-log(p)_comb");
      if (qvalOpt)
        gzprintf(log.gzf, "\t-log(q)");
      gzprintf(log.gzf, "\tsignif\n");
    } else {
      fprintf(log.f, "chr\tstart\tend");
      for (int i = 0; i < n; i++)
        fprintf(log.f, "\t-log(p)_%d", i);
      fprintf(log.f, "\t-log(p)_comb");
      if (qvalOpt)
        fprintf(log.f, "\t-log(q)");
      fprintf(log.f, "\tsignif\n");
    }
  } else {
    // single sample: logfile has pileups and p-/q-values
    if (gzOut) {
      gzprintf(log.gzf, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)");
      if (qvalOpt)
        gzprintf(log.gzf, "\t-log(q)");
      gzprintf(log.gzf, "\tsignif\n");
    } else {
      fprintf(log.f, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)");
      if (qvalOpt)
        fprintf(log.f, "\t-log(q)");
      fprintf(log.f, "\tsignif\n");
    }
  }
}

/* void printIntervalN()
 * Print bedgraph(ish) interval for multiple replicates.
 *   Values: -log(p) for each replicate, combined -log(p),
 *   -log(q), and significance ('*') for each.
 */
void printIntervalN(File out, bool gzOut, char* name,
    uint32_t start, uint32_t end, Pileup** p, int n,
    uint32_t idx[], float pval, float qval, bool sig) {
  if (gzOut) {
    gzprintf(out.gzf, "%s\t%d\t%d", name, start, end);
    for (int i = 0; i < n; i++)
      if (p[i] == NULL)
        gzprintf(out.gzf, "\tn/a");
      else
        gzprintf(out.gzf, "\t%f", p[i]->cov[idx[i]]);
    gzprintf(out.gzf, "\t%f", pval);
    if (qval != -1.0f)
      gzprintf(out.gzf, "\t%f", qval);
    gzprintf(out.gzf, "%s\n", sig ? "\t*" : "");
  } else {
    fprintf(out.f, "%s\t%d\t%d", name, start, end);
    for (int i = 0; i < n; i++)
      if (p[i] == NULL)
        fprintf(out.f, "\tn/a");
      else
        fprintf(out.f, "\t%f", p[i]->cov[idx[i]]);
    fprintf(out.f, "\t%f", pval);
    if (qval != -1.0f)
      fprintf(out.f, "\t%f", qval);
    fprintf(out.f, "%s\n", sig ? "\t*" : "");
  }
}

/* void printInterval()
 * Print bedgraph(ish) interval for a single replicate.
 *   Values: pileups (treatment and control), -log(p),
 *   -log(q), and significance ('*') for each.
 */
void printInterval(File out, bool gzOut, char* name,
    uint32_t start, uint32_t end, float treatVal,
    float ctrlVal, float pval, float qval, bool sig) {
  if (gzOut) {
    gzprintf(out.gzf, "%s\t%d\t%d\t%f\t%f\t%f",
      name, start, end, treatVal, ctrlVal, pval);
    if (qval != -1.0f)
      gzprintf(out.gzf, "\t%f", qval);
    gzprintf(out.gzf, "%s\n", sig ? "\t*" : "");
  } else {
    fprintf(out.f, "%s\t%d\t%d\t%f\t%f\t%f",
      name, start, end, treatVal, ctrlVal, pval);
    if (qval != -1.0f)
      fprintf(out.f, "\t%f", qval);
    fprintf(out.f, "%s\n", sig ? "\t*" : "");
  }
}

/* void printLog()
 * Control printing of stats for an interval.
 */
void printLog(File log, bool gzOut, Chrom* chr,
    uint32_t start, int n, uint32_t m, uint32_t j,
    uint32_t k, uint32_t idx[], bool qvalOpt,
    bool sig) {
  if (! n) {
    // single replicate
    printInterval(log, gzOut, chr->name,
      start, chr->pval[n]->end[m],
      chr->treat->cov[j], chr->ctrl->cov[k],
      chr->pval[n]->cov[m],
      qvalOpt ? chr->qval->cov[m] : -1.0f, sig);
  } else {
    // multiple replicates
    printIntervalN(log, gzOut, chr->name,
      start, chr->pval[n]->end[m],
      chr->pval, n, idx, chr->pval[n]->cov[m],
      qvalOpt ? chr->qval->cov[m] : -1.0f, sig);
    // update indexes into pval arrays
    for (int r = 0; r < n; r++)
      if (chr->pval[r] != NULL
          && chr->pval[r]->end[idx[r]] == chr->pval[n]->end[m])
        idx[r]++;
  }
}

/* void printPeak()
 * Print peaks in ENCODE narrowPeak format.
 */
void printPeak(File out, bool gzOut, char* name,
    int64_t start, int64_t end, int count, float val,
    float fe, float pval, float qval, uint32_t pos) {
  if (gzOut)
    gzprintf(out.gzf, "%s\t%ld\t%ld\tpeak_%d\t%d\t.\t%f\t%f\t%f\t%d\n",
      name, start, end, count,
      MIN((unsigned int) (val * 10.0f + 0.5f), 1000),
      fe, pval, qval, pos);
  else
    fprintf(out.f, "%s\t%ld\t%ld\tpeak_%d\t%d\t.\t%f\t%f\t%f\t%d\n",
      name, start, end, count,
      MIN((unsigned int) (val * 10.0f + 0.5f), 1000),
      fe, pval, qval, pos);
}

/* int callPeaks()
 * Call peaks. Produce output on the fly.
 *   Return number of peaks.
 */
int callPeaks(File out, File log, bool logOpt, bool gzOut,
    Chrom* chrom, int chromLen, int n, float pqvalue,
    bool qvalOpt, int minLen, int maxGap) {

  if (logOpt)
    printLogHeader(log, gzOut, n, qvalOpt);

  // loop through chroms
  int count = 0;      // count of peaks
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || (qvalOpt && chr->qval == NULL)
        || (! qvalOpt && chr->pval[n] == NULL) )
      continue;

    // create indexes into arrays for logging purposes
    //   (treat/ctrl pileup [if n == 0] and p-value arrays [if n > 0])
    uint32_t j = 0, k = 0;  // indexes into chr->treat, chr->ctrl
    uint32_t idx[n];        // indexes into each pval array
    for (int r = 0; r < n; r++)
      idx[r] = 0;

    // reset peak variables
    int64_t peakStart = -1, peakEnd = -1; // ends of potential peak
    float summitVal = -1.0f;              // summit p/q value
    uint32_t summitPos = 0;               // distance from peakStart to summit
    uint32_t summitLen = 0;               // length of summit interval
    float summitPval = -1.0f, summitQval = -1.0f; // summit p- and q-values
    float summitFE = -1.0f;               // summit fold enrichment

    // loop through intervals (defined by chr->pval[n])
    uint32_t start = 0;    // start of interval
    for (uint32_t m = 0; m < chr->pvalLen[n]; m++) {

      bool sig = false;
      float val = qvalOpt ? chr->qval->cov[m] : chr->pval[n]->cov[m];
      if ( val > pqvalue ) {

        // interval reaches significance
        sig = true;
        if (peakStart == -1)
          peakStart = start;  // start new potential peak
        peakEnd = chr->pval[n]->end[m];  // end of potential peak

        // check if interval is summit for this peak
        if (val > summitVal) {
          summitVal = val;
          if (! n)
            summitFE = chr->ctrl->cov[k] ?
              chr->treat->cov[j] / chr->ctrl->cov[k] : FLT_MAX;
          summitPval = chr->pval[n]->cov[m];
          summitQval = qvalOpt ? chr->qval->cov[m] : -1.0f;
          summitPos = (peakEnd + (m ? chr->pval[n]->end[m-1] : 0) ) / 2
            - peakStart;  // midpoint of interval
          summitLen = peakEnd - (m ? chr->pval[n]->end[m-1] : 0);
        } else if (val == summitVal) {
          // update summitPos only if interval is longer
          uint32_t len = chr->pval[n]->end[m]
            - (m ? chr->pval[n]->end[m-1] : 0);
          if (len > summitLen) {
            summitPos = (peakEnd + (m ? chr->pval[n]->end[m-1] : 0) ) / 2
              - peakStart;  // midpoint of interval
            summitLen = len;
          }
        }

      } else {

        // interval does not reach significance
        if (peakStart != -1
            && chr->pval[n]->end[m] - peakEnd > maxGap) {
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

      // print stats for interval
      if (logOpt)
        printLog(log, gzOut, chr, start, n, m, j, k,
          idx, qvalOpt, sig);

      // update chr->treat and chr->ctrl indexes
      if (! n) {
        if (chr->ctrl->end[k] < chr->treat->end[j])
          k++;
        else {
          if (chr->ctrl->end[k] == chr->treat->end[j])
            k++;
          j++;
        }
      }

      start = chr->pval[n]->end[m];
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
    Chrom* chrom, int chromLen, int* sample, float pqvalue,
    bool qvalOpt, int minLen, int maxGap, bool verbose) {

  // calculate combined p-values for multiple replicates
  if (*sample > 1) {
    combinePval(chrom, chromLen, *sample);
    (*sample)++;
  }

  // calculate genome length (only chroms that are not
  //   skipped and have had p-values calculated)
  uint64_t genomeLen = 0;
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip && chr->pval[*sample - 1] != NULL) {
      genomeLen += chr->len;
      for (int j = 0; j < chr->bedLen; j++)
        genomeLen -= chr->bedEnd[j] - chr->bedSt[j];
    }
  }

  if (verbose) {
    fprintf(stderr, "Peak-calling parameters:\n");
    fprintf(stderr, "  Genome length: %ldbp\n", genomeLen);
    fprintf(stderr, "  Significance threshold: -log(%c) > %.3f\n",
      qvalOpt ? 'q' : 'p', pqvalue);
    fprintf(stderr, "  Max. gap between sites: %dbp\n", maxGap);
    fprintf(stderr, "  Min. peak length: %dbp\n", minLen);
  }

  // compute q-values
  if (qvalOpt)
    computeQval(chrom, chromLen, genomeLen, *sample - 1);

  // identify peaks
  int count = callPeaks(out, log, logOpt, gzOut, chrom,
    chromLen, *sample - 1, pqvalue, qvalOpt, minLen,
    maxGap);
  if (verbose)
    fprintf(stderr, "Peaks identified: %d\n", count);
}

/*** Calculate p-values ***/

/* uint32_t countIntervals()
 * Count the number of pileup intervals to create
 *   for a composite.
 */
uint32_t countIntervals(Chrom* chr) {
  uint32_t num = 0;
  uint32_t k = 0;
  for (uint32_t j = 0; j < chr->treatLen; j++) {
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

/* void printPileHeader()
 * Print header of the bedgraph-ish pileup log file.
 */
void printPileHeader(File pile, char* treatName,
    char* ctrlName, bool gzOut) {
  if (gzOut) {
    gzprintf(pile.gzf, "# treatment file: %s; control file: %s\n",
      treatName, ctrlName && strcmp(ctrlName, "null") ? ctrlName : "n/a");
    gzprintf(pile.gzf, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)\n");
  } else {
    fprintf(pile.f, "# treatment file: %s; control file: %s\n",
      treatName, ctrlName && strcmp(ctrlName, "null") ? ctrlName : "n/a");
    fprintf(pile.f, "chr\tstart\tend\ttreatment\tcontrol\t-log(p)\n");
  }
}

/* void printPile()
 * Print bedgraph-ish interval of treatment/control
 *   pileup values and p-value.
 */
void printPile(File pile, char* name, uint32_t start,
    uint32_t end, float treat, float ctrl, float pval,
    bool gzOut) {
  if (gzOut)
    gzprintf(pile.gzf, "%s\t%d\t%d\t%f\t%f\t%f\n",
      name, start, end, treat, ctrl, pval);
  else
    fprintf(pile.f, "%s\t%d\t%d\t%f\t%f\t%f\n",
      name, start, end, treat, ctrl, pval);
}

/* float calcPval()
 * Calculate -log10(p) using an exponential distribution
 *   with parameter beta (ctrlVal) and observation treatVal.
 */
float calcPval(float treatVal, float ctrlVal) {
  if (ctrlVal == 0.0f)
    return treatVal == 0.0f ? 0.0f : FLT_MAX;
  double pval = treatVal / ctrlVal * M_LOG10E;
  return pval > FLT_MAX ? FLT_MAX : (float) pval;
}

/* void savePval()
 * Create and save p-values as pileups for each chrom.
 */
void savePval(Chrom* chrom, int chromLen, int n,
    File pile, bool pileOpt, char* treatName,
    char* ctrlName, bool gzOut) {

  // print log header
  if (pileOpt)
    printPileHeader(pile, treatName, ctrlName, gzOut);

  // create pileups for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip)
      continue;

    // fill in missing pval arrays from previous samples
    if (chr->sample < n) {
      chr->pval = (Pileup**) memrealloc(chr->pval,
        n * sizeof(Pileup*));
      for (int j = n - 1; j >= chr->sample; j--)
        chr->pval[j] = NULL;
      chr->sample = n;
    }

    // p-values not to be saved: append a NULL
    if (! chr->save) {
      chr->pval = (Pileup**) memrealloc(chr->pval,
        (n + 1) * sizeof(Pileup*));
      chr->pval[n] = NULL;
      chr->sample++;
      continue;
    }

    // create 'pileup' arrays for p-values
    uint32_t num = countIntervals(chr);
    chr->pval = (Pileup**) memrealloc(chr->pval,
      (n + 1) * sizeof(Pileup*));
    chr->pval[n] = (Pileup*) memalloc(sizeof(Pileup));
    chr->pval[n]->end = (uint32_t*) memalloc(num * sizeof(uint32_t));
    chr->pval[n]->cov = (float*) memalloc(num * sizeof(float));
    chr->pvalLen = (uint32_t*) memrealloc(chr->pvalLen,
      (n + 1) * sizeof(uint32_t));
    chr->pvalLen[n] = num;
    chr->sample++;

    // save p-values to arrays
    Pileup* p = chr->pval[n];
    uint32_t start = 0;    // start of interval
    uint32_t j = 0, k = 0;
    for (uint32_t m = 0; m < num; m++) {
      if (chr->ctrl->end[k] < chr->treat->end[j]) {
        p->end[m] = chr->ctrl->end[k];
        p->cov[m] = calcPval(chr->treat->cov[j],
          chr->ctrl->cov[k]);
        if (pileOpt)
          printPile(pile, chr->name, start, p->end[m],
            chr->treat->cov[j], chr->ctrl->cov[k],
            p->cov[m], gzOut);
        k++;
      } else {
        p->end[m] = chr->treat->end[j];
        p->cov[m] = calcPval(chr->treat->cov[j],
          chr->ctrl->cov[k]);
        if (pileOpt)
          printPile(pile, chr->name, start, p->end[m],
            chr->treat->cov[j], chr->ctrl->cov[k],
            p->cov[m], gzOut);
        if (chr->ctrl->end[k] == chr->treat->end[j])
          k++;
        j++;
      }
      start = p->end[m];
    }

  }
}

/*** Save treatment/control pileup values ***/

/* void saveConst()
 * Save a given value as pileup for a full chromosome.
 */
void saveConst(Pileup* p, uint32_t* size, uint32_t* mem,
    uint32_t len, float val) {
  if (! *mem) {
    p->end = (uint32_t*) memalloc(sizeof(uint32_t));
    p->cov = (float*) memalloc(sizeof(float));
    *mem = 1;
  }
  p->end[0] = len;
  p->cov[0] = val;
  *size = 1;
}

/* float calcLambda()
 * Calculate a background lambda value: sum of fragment
 *   lengths divided by total genome length.
 */
float calcLambda(Chrom* chrom, int chromLen,
    double fragLen, uint64_t* genomeLen) {
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip && chr->save) {
      *genomeLen += chr->len;
      for (int j = 0; j < chr->bedLen; j++)
        *genomeLen -= chr->bedEnd[j] - chr->bedSt[j];
    }
  }
  if (! *genomeLen)
    exit(error("", ERRGEN));
  return fragLen / *genomeLen;
}

/* void savePileupNoCtrl()
 * When no control is available, save the control
 *   pileup as the background lambda value.
 */
void savePileupNoCtrl(Chrom* chrom, int chromLen,
    double fragLen, uint64_t* genomeLen) {
  float lambda = calcLambda(chrom, chromLen, fragLen,
    genomeLen);
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save)
      continue;
    saveConst(chr->ctrl, &chr->ctrlLen, &chr->ctrlMem,
      chr->len, lambda);
  }
}

/* float getVal()
 * Reconstruct a float value from the given
 *   int and 8-bit encoded fractional part.
 */
float getVal(int32_t cov, uint8_t frac) {
  return (float) cov
    + ((frac & 0x7) / 8.0f)
    + (((frac >> 3) & 0x3) / 6.0f)
    + (((frac >> 5) & 0x7) / 10.0f);
}

/* float updateVal()
 * Update pileup value (int cov and 8-bit encoded
 *   fractional part) by adding the given int (dCov)
 *   and fractional part (dFrac) from the 'diff' arrays.
 *   Return reconstructed value via getVal().
 */
float updateVal(int16_t dCov, uint8_t dFrac, int32_t* cov,
    uint8_t* frac) {

  // add ints
  *cov += dCov;
  if (! dFrac) {
//fprintf(stderr, "  adding %d     -> %d / %02x", dCov, *cov, *frac);
//while (!getchar()) ;
    if (*cov < 0)
      exit(error(errMsg[ERRPILE], ERRISSUE));
    return getVal(*cov, *frac);
  }

  // sum eighths (and collect halves)
  int half = 0;
  if (dFrac & 0x7) {
    int sum8 = (dFrac & 0x7) + (*frac & 0x7);
    while (sum8 > 3) {
      half++;
      sum8 -= 4;
    }
    *frac = (*frac & 0xF8) | sum8;
  } else if (*frac & 0x4) {
    half++;
    *frac &= 0xFB;
  }

  // sum sixths
  if (dFrac & 0x18) {
    int sum6 = ((dFrac >> 3) & 0x3) + ((*frac >> 3) & 0x3);
    if (sum6 > 2) {
      half++;
      sum6 -= 3;
    }
    *frac = (*frac & 0xE7) | (sum6 << 3);
  }

  // sum tenths
  if (dFrac & 0xE0) {
    int sum10 = ((dFrac >> 5) & 0x7) + ((*frac >> 5) & 0x7);
    if (sum10 > 4) {
      half++;
      sum10 -= 5;
    }
    *frac = (*frac & 0x1F) | (sum10 << 5);
  }

  // combine halves
  while (half > 1) {
    (*cov)++;
    half -= 2;
  }
  if (half)
    *frac |= 0x4;

  // check for negative pileup
  if (*cov < 0)
    exit(error(errMsg[ERRPILE], ERRISSUE));

//fprintf(stderr, "  adding %d / %02x -> %d / %02x", dCov, dFrac, *cov, *frac);
//while (!getchar()) ;
  return getVal(*cov, *frac);
}

/* float calcFactor
 * Calculate the scaling factor of treatment fragment
 *   lengths to control. Also set ctrlLen for each
 *   Chrom* (to be corrected in savePileupCtrl()).
 */
float calcFactor(Chrom* chrom, int chromLen,
    double fragLen) {

  // sum weighted lengths of control fragments
  double ctrlFrag = 0.0;
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save || chr->diff == NULL)
      continue;

    // initialize variables
    uint32_t num = 1;       // number of intervals to create for pileup
    Diff* d = chr->diff;
    int32_t cov = 0;        // current pileup value
    uint8_t frac = 0;       // current pileup value (fraction part)
    float val = updateVal(d->cov[0], d->frac[0], &cov, &frac);

    // calculate fragment lengths along the chrom
    uint32_t start = 0;     // beginning coordinate of interval
    uint32_t j;
    for (j = 1; j < chr->len; j++) {
      if (d->cov[j] || d->frac[j]) {

        // sum fragment length (weighted by val)
        ctrlFrag += (j - start) * val;
        start = j;

        // update pileup value
        val = updateVal(d->cov[j], d->frac[j], &cov, &frac);
        num++;

      }
    }

    // save final interval
    ctrlFrag += (j - start) * val;
    chr->ctrlLen = num; // save number of intervals
  }

  // return ratio of treatment frags to ctrl frags
  if (! ctrlFrag)
    return 1.0f;
  return fragLen / ctrlFrag;
}

/* void savePileupCtrl()
 * Save pileup values for control sample(s) from
 *   'diff' arrays and background lambda value.
 */
void savePileupCtrl(Chrom* chrom, int chromLen,
    double fragLen, uint64_t* genomeLen, bool verbose) {

  // calculate background lambda value
  float lambda = calcLambda(chrom, chromLen, fragLen,
    genomeLen);

  // calculate scale factor (treatment / control)
  float factor = calcFactor(chrom, chromLen, fragLen);
  if (verbose) {
    fprintf(stderr, "  Scaling factor for control pileup: %f\n", factor);
    if (factor > 5.0f)
      fprintf(stderr, "  ** Warning! Large scaling may mask true signal **\n");
  }

  // create pileup for each chrom
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save)
      continue;

    // if no read coverage, save constant pileup of lambda
    if (chr->diff == NULL) {
      saveConst(chr->ctrl, &chr->ctrlLen, &chr->ctrlMem,
        chr->len, lambda);
      continue;
    }

    // expand pileup arrays (if necessary)
    if (chr->ctrlLen > chr->ctrlMem) {
      chr->ctrl->end = (uint32_t*) memrealloc(chr->ctrl->end,
        chr->ctrlLen * sizeof(uint32_t));
      chr->ctrl->cov = (float*) memrealloc(chr->ctrl->cov,
        chr->ctrlLen * sizeof(float));
      chr->ctrlMem = chr->ctrlLen;
    }

    // initialize pileup values
    Diff* d = chr->diff;
    int32_t cov = 0;      // current pileup value
    uint8_t frac = 0;     // current pileup value (fraction part)
    float val = factor * updateVal(d->cov[0], d->frac[0],
      &cov, &frac);
    float net = MAX(val, lambda);

    // save pileup values along the chrom
    uint32_t pos = 0;     // position in pileup arrays
    uint32_t j;
    for (j = 1; j < chr->len; j++) {
      if (d->cov[j] || d->frac[j]) {

        // update pileup value
        val = factor * updateVal(d->cov[j], d->frac[j],
          &cov, &frac);

        // save interval if value has changed
        if (net != MAX(val, lambda)) {
          chr->ctrl->end[pos] = j;
          chr->ctrl->cov[pos] = net;
          pos++;
          net = MAX(val, lambda);
        }

      }
    }

    // save final interval
    chr->ctrl->end[pos] = j;
    chr->ctrl->cov[pos] = net;

    // update array length
    if (pos >= chr->ctrlLen) {
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "%s (%s)", errMsg[ERRARRC], chr->name);
      exit(error(msg, ERRISSUE));
    }
    chr->ctrlLen = pos + 1;

    // check for val error (should end at 0)
    val = updateVal(d->cov[j], d->frac[j], &cov, &frac);
    if (val) {
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "Control pileup for ref %s finishes at %f (not 0.0)",
        chr->name, val);
      exit(error(msg, ERRISSUE));
    }
  }

}

/* double savePileupTreat()
 * Save pileup values for treatment sample(s) from
 *   'diff' arrays.
 *   Return total length of all fragments (weighted).
 */
double savePileupTreat(Chrom* chrom, int chromLen) {

  // create pileup for each chrom
  double fragLen = 0.0;  // weighted fragment length
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save)
      continue;

    // if no read coverage, save constant pileup of 0
    if (chr->diff == NULL) {
      saveConst(chr->treat, &chr->treatLen, &chr->treatMem,
        chr->len, 0.0f);
      continue;
    }

//fprintf(stderr, "chrom %s\n", chr->name);

    // determine number of pileup intervals
    Diff* d = chr->diff;
    uint32_t num = 1;
    for (uint32_t j = 1; j < chr->len; j++)
      if (d->cov[j] || d->frac[j])
        num++;

    // expand pileup arrays (if necessary)
    if (num > chr->treatMem) {
      chr->treat->end = (uint32_t*) memrealloc(chr->treat->end,
        num * sizeof(uint32_t));
      chr->treat->cov = (float*) memrealloc(chr->treat->cov,
        num * sizeof(float));
      chr->treatMem = num;
    }
    chr->treatLen = num;

    // initialize pileup values
    int32_t cov = 0;      // current pileup value
    uint8_t frac = 0;     // current pileup value (fraction part)
    float val = updateVal(d->cov[0], d->frac[0], &cov, &frac);

    // save pileup values along the chrom
    uint32_t start = 0;   // beginning coordinate of interval
    uint32_t pos = 0;     // position in pileup arrays
    uint32_t j;
    for (j = 1; j < chr->len; j++) {
      if (d->cov[j] || d->frac[j]) {

        // save end of interval and pileup value
        chr->treat->end[pos] = j;
        chr->treat->cov[pos] = val;
        pos++;

        // keep track of fragment length (weighted by val)
        fragLen += (j - start) * val;
        start = j;

        // update pileup value
        val = updateVal(d->cov[j], d->frac[j], &cov, &frac);

//fprintf(stderr, "  pos %d: %.9f -> %.9f", j, getVal(d->cov[j], d->frac[j]), val);
//while(!getchar()) ;
      }
    }

    // save final interval
    chr->treat->end[pos] = j;
    chr->treat->cov[pos] = val;
    fragLen += (j - start) * val;

    // verify array length
    if (pos + 1 != chr->treatLen) {
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "%s (%s)", errMsg[ERRARR], chr->name);
      exit(error(msg, ERRISSUE));
    }

    // check for val error (should end at 0)
    val = updateVal(d->cov[j], d->frac[j], &cov, &frac);
    if (val) {
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "Treatment pileup for ref %s finishes at %f (not 0.0)",
        chr->name, val);
      exit(error(msg, ERRISSUE));
    }
  }

  if (fragLen == 0.0)
    exit(error("", ERRTREAT));
  return fragLen;
}

/*** Fractional alignment accounting ***/

/* void addFrac()
 * Add a fractional count (1/count) to *frac,
 *   which has this 8-bit encoding:
 *       000     00      000
 *     tenths  sixths  eighths
 *   Carry over is applied to *cov.
 *
 * This function (as well as subFrac(), below) is
 *   written at the bit level to optimize efficiency,
 *   and hence is nearly impossible to debug/maintain.
 *   Apologies in advance to those attempting to do so.
 */
void addFrac(int16_t* cov, uint8_t* frac, uint8_t count) {
  switch (count) {
    case 8:
      if ((*frac & 0x7) == 0x7) {
        (*cov)++;
        *frac &= 0xF8;
      } else
        (*frac)++;
      break;
    case 4:
      if ((*frac & 0x6) == 0x6) {
        (*cov)++;
        *frac &= 0xF9;
      } else
        *frac += 0x2;
      break;
    case 2:
      if (*frac & 0x4) {
        (*cov)++;
        *frac &= 0xFB;
      } else
        *frac |= 0x4;
      break;
    case 6:
      if (*frac & 0x10) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0xEB;
        } else {
          *frac |= 0x4;
          *frac &= 0xEF;
        }
      } else
        *frac += 0x8;
      break;
    case 3:
      if (*frac & 0x8) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0xF3;
        } else {
          *frac |= 0x4;
          *frac &= 0xF7;
        }
      } else if (*frac & 0x10) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0xEB;
        } else {
          *frac |= 0x4;
          *frac &= 0xEF;
        }
        *frac |= 0x8;
      } else
        *frac |= 0x10;
      break;
    case 10:
      if (*frac & 0x80) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0x7B;
        } else {
          *frac |= 0x4;
          *frac &= 0x7F;
        }
      } else
        *frac += 0x20;
      break;
    case 5:
      if (*frac & 0x80) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0x7B;
        } else {
          *frac |= 0x4;
          *frac &= 0x7F;
        }
        *frac += 0x20;
      } else if ((*frac & 0x60) == 0x60) {
        if (*frac & 0x4) {
          (*cov)++;
          *frac &= 0x9B;
        } else {
          *frac |= 0x4;
          *frac &= 0x9F;
        }
      } else
        *frac += 0x40;
      break;
    default: ;
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "%s (%d)", errMsg[ERRALNS], count);
      exit(error(msg, ERRISSUE));
  }
}

/* void subFrac()
 * Subtract a fractional count to *frac.
 *   See addFrac() above for a description of
 *   the 8-bit encoding.
 */
void subFrac(int16_t* cov, uint8_t* frac, uint8_t count) {
  switch (count) {
    case 8:
      if (*frac & 0x7)
        (*frac)--;
      else {
        (*cov)--;
        *frac |= 0x7;
      }
      break;
    case 4:
      if (*frac & 0x6)
        *frac -= 0x2;
      else {
        (*cov)--;
        *frac |= 0x6;
      }
      break;
    case 2:
      if (*frac & 0x4)
        *frac &= 0xFB;
      else {
        (*cov)--;
        *frac |= 0x4;
      }
      break;
    case 6:
      if (*frac & 0x18)
        *frac -= 0x8;
      else if (*frac & 0x4) {
        *frac |= 0x10;
        *frac &= 0xFB;
      } else {
        (*cov)--;
        *frac |= 0x14;
      }
      break;
    case 3:
      if (*frac & 0x10)
        *frac &= 0xEF;
      else if (*frac & 0x4) {
        *frac &= 0xFB;
        *frac += 0x8;
      } else {
        (*cov)--;
        *frac += 0xC;
      }
      break;
    case 10:
      if (*frac & 0xE0)
        *frac -= 0x20;
      else if (*frac & 0x4) {
        *frac &= 0xFB;
        *frac |= 0x80;
      } else {
        (*cov)--;
        *frac |= 0x84;
      }
      break;
    case 5:
      if (*frac & 0xC0)
        *frac -= 0x40;
      else if (*frac & 0x4) {
        *frac &= 0xFB;
        *frac += 0x60;
      } else {
        (*cov)--;
        *frac |= 0x4;
        *frac += 0x60;
      }
      break;
    default: ;
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "%s (%d)", errMsg[ERRALNS], count);
      exit(error(msg, ERRISSUE));
  }
}

/*** Convert alignments to intervals ***/

/* void printBED()
 * Print a BED interval for a read/fragment.
 *   Append the aln count, 'C'ontrol/'T'reatment, and
 *   sample number to the read name (4th column).
 */
void printBED(File bed, bool gzOut, char* chr,
    int64_t start, int64_t end, char* qname,
    uint8_t count, bool ctrl, int sample) {
  if (gzOut)
    gzprintf(bed.gzf, "%s\t%ld\t%ld\t%s_%d_%c_%d\n",
      chr, start, end, qname, count, ctrl ? 'C' : 'T',
      sample);
  else
    fprintf(bed.f, "%s\t%ld\t%ld\t%s_%d_%c_%d\n",
      chr, start, end, qname, count, ctrl ? 'C' : 'T',
      sample);
}

/* uint32_t saveInterval()
 * Check the validity of the start/end coordinates
 *   of a read/fragment. Save the ends to the 'diff'
 *   arrays of the given Chrom*.
 *   Return the fragment length.
 */
uint32_t saveInterval(Chrom* c, int64_t start, int64_t end,
    char* qname, uint8_t count, File bed, bool bedOpt,
    bool gzOut, bool ctrl, int sample, bool verbose) {

  // check validity of positions
  if (start < 0) {
    if (verbose)
      fprintf(stderr, "Warning! Read %s prevented from extending below 0 on %s\n",
        qname, c->name);
    start = 0;
  }
  if (start >= c->len) {
    char* msg = (char*) memalloc(MAX_ALNS);
    sprintf(msg, "Read %s, ref. %s", qname, c->name);
    exit(error(msg, ERRPOS));
  }
  if (end > c->len) {
    if (verbose)
      fprintf(stderr, "Warning! Read %s prevented from extending past %d on %s\n",
        qname, c->len, c->name);
    end = c->len;
  }

  // create 'diff' arrays if necessary
  if (c->diff == NULL) {
    c->diff = (Diff*) memalloc(sizeof(Diff));
    c->diff->frac = (uint8_t*) memalloc((1 + c->len) * sizeof(uint8_t));
    c->diff->cov = (int16_t*) memalloc((1 + c->len) * sizeof(int16_t));
    for (int i = 0; i < 1 + c->len; i++) {
      c->diff->frac[i] = 0;
      c->diff->cov[i] = 0;
    }
  }

  // check for overflow/underflow (c->diff->cov is int16_t)
  if (c->diff->cov[start] == SHRT_MAX) {
    if (verbose)
      fprintf(stderr, "Warning! Read %s, alignment at (%s, %ld-%ld) skipped due to overflow\n",
        qname, c->name, start, end);
    return 0;
  }
  if (c->diff->cov[end] == SHRT_MIN) {
    if (verbose)
      fprintf(stderr, "Warning! Read %s, alignment at (%s, %ld-%ld) skipped due to underflow\n",
        qname, c->name, start, end);
    return 0;
  }

  // add counts to diff array(s)
  if (count == 1) {
    c->diff->cov[start]++;
    c->diff->cov[end]--;
  } else {
    // add fractional count
    addFrac(&c->diff->cov[start], &c->diff->frac[start], count);
    subFrac(&c->diff->cov[end], &c->diff->frac[end], count);
  }

  // print BED interval
  if (bedOpt)
    printBED(bed, gzOut, c->name, start, end, qname,
      count, ctrl, sample);

  return end - start;
}

/* void processAvgExt()
 * Save complete intervals for unpaired alignments
 *   with "extend to average length" option, after
 *   calculating average length from paired alns.
 */
void processAvgExt(Aln** unpair, int unpairIdx,
    int unpairLen, double totalLen, int pairedPr,
    File bed, bool bedOpt, bool gzOut, bool ctrl,
    int sample, bool verbose) {

  // determine average fragment length
  int avgLen = 0;
  if (! pairedPr && verbose) {
    fprintf(stderr, "Warning! No paired alignments to calculate avg ");
    fprintf(stderr, "frag length --\n  Printing singletons \"as is\"\n");
  } else
    avgLen = (int) (totalLen / pairedPr + 0.5);

  // process each alignment
  for (int i = 0; i <= unpairIdx; i++) {

    int end = (i == unpairIdx ? unpairLen : MAX_SIZE);
    for (int j = 0; j < end; j++) {

      Aln* a = unpair[i] + j;
      if (! avgLen)
        saveInterval(a->chrom, a->pos[0], a->pos[1], a->name,
          a->count, bed, bedOpt, gzOut, ctrl, sample,
          verbose);
      else if (a->strand)
        saveInterval(a->chrom, a->pos[0], a->pos[0] + avgLen,
          a->name, a->count, bed, bedOpt, gzOut, ctrl,
          sample, verbose);
      else
        saveInterval(a->chrom, (signed) (a->pos[1] - avgLen),
          a->pos[1], a->name, a->count, bed, bedOpt, gzOut,
          ctrl, sample, verbose);

      // free memory
      free(a->name);
    }

  }
}

/* void saveAvgExt()
 * Save info for an unpaired alignment to list
 *   (for "extend to average length" option), for
 *   later processing by saveAvgExt().
 */
void saveAvgExt(char* qname, Aln* b, uint8_t count,
    Aln*** unpair, int* unpairIdx, int* unpairLen,
    int* unpairMem) {

  // alloc memory if necessary
  if (*unpairLen == 0 && *unpairIdx == *unpairMem) {
    *unpair = (Aln**) memrealloc(*unpair,
      (*unpairMem + 1) * sizeof(Aln*));
    (*unpair)[*unpairMem] = (Aln*) memalloc(MAX_SIZE
      * sizeof(Aln));
    (*unpairMem)++;
  }

  // copy alignment info
  Aln* a = (*unpair)[*unpairIdx] + *unpairLen;
  a->name = (char*) memalloc(1 + strlen(qname));
  strcpy(a->name, qname);
  a->chrom = b->chrom;
  a->strand = b->strand;
  a->pos[0] = b->pos[0];
  a->pos[1] = b->pos[1];
  a->count = count;

  (*unpairLen)++;
  if (*unpairLen == MAX_SIZE) {
    *unpairLen = 0;
    (*unpairIdx)++;
  }
}

/* void saveSingle()
 * Control processing of singleton alignments
 *   (either keeping them as is, or extending
 *   to a given length).
 */
void saveSingle(char* qname, Aln* a, uint8_t count,
    bool extendOpt, int extend, bool atacOpt,
    int atacLen5, int atacLen3, File bed, bool bedOpt,
    bool gzOut, bool ctrl, int sample, bool verbose) {
  if (extendOpt) {
    if (a->strand)
      saveInterval(a->chrom, a->pos[0], a->pos[0] + extend,
        qname, count, bed, bedOpt, gzOut, ctrl, sample,
        verbose);
    else
      saveInterval(a->chrom, (signed) (a->pos[1] - extend),
        a->pos[1], qname, count, bed, bedOpt, gzOut, ctrl,
        sample, verbose);
  } else if (atacOpt) {
    if (a->strand)
      saveInterval(a->chrom, (signed) (a->pos[0] - atacLen5),
        a->pos[0] + atacLen3, qname, count, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
    else
      saveInterval(a->chrom, (signed) (a->pos[1] - atacLen3),
        a->pos[1] + atacLen5, qname, count, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
  } else
    saveInterval(a->chrom, a->pos[0], a->pos[1], qname,
      count, bed, bedOpt, gzOut, ctrl, sample, verbose);
}

/* uint32_t saveFragAtac()
 * In ATAC-seq mode, save intervals for each end of a full
 *   fragment. If they overlap, just save one big interval.
 *   Return total length.
 */
uint32_t saveFragAtac(Chrom* c, uint32_t start,
    uint32_t end, int atacLen5, int atacLen3,
    char* qname, uint8_t count, File bed, bool bedOpt,
    bool gzOut, bool ctrl, int sample, bool verbose) {
  if (start + atacLen3 >= (signed) (end - atacLen3))
    // expanded intervals overlap: just save one
    return saveInterval(c, (signed) (start - atacLen5),
      end + atacLen5, qname, count, bed, bedOpt, gzOut,
      ctrl, sample, verbose);
  // save two intervals
  return saveInterval(c, (signed) (start - atacLen5),
      start + atacLen3, qname, count, bed, bedOpt,
      gzOut, ctrl, sample, verbose)
    + saveInterval(c, (signed) (end - atacLen3),
      end + atacLen5, qname, count, bed, bedOpt,
      gzOut, ctrl, sample, verbose);
}

/* uint32_t saveFragment()
 * Save full fragment for a proper pair. Return length.
 */
uint32_t saveFragment(char* qname, Aln* a, uint8_t count,
    bool atacOpt, int atacLen5, int atacLen3,
    File bed, bool bedOpt, bool gzOut, bool ctrl,
    int sample, bool verbose) {
  // ensure start < end
  uint32_t start, end;
  if (a->pos[0] > a->pos[1]) {
    start = a->pos[1];
    end = a->pos[0];
  } else {
    start = a->pos[0];
    end = a->pos[1];
  }
  if (atacOpt)
    return saveFragAtac(a->chrom, start, end, atacLen5,
      atacLen3, qname, count, bed, bedOpt, gzOut, ctrl,
      sample, verbose);
  return saveInterval(a->chrom, start, end, qname,
    count, bed, bedOpt, gzOut, ctrl, sample, verbose);
}

/*** Process a set of alignments ***/

/* void subsampleSingle()
 * For sets of single alns at an invalid count (>10, 9, 7),
 *   find a more stringent score.
 */
void subsampleSingle(Aln* aln, int alnLen, bool first,
    uint8_t* count, float* score) {

  // insertion sort of aln scores
  float arr[*count];  // sorted array of aln scores
  int k = 0;          // count of valid alns analyzed
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (! a->paired && a->first == first
        && a->score >= *score && a->chrom->save
        && ! a->chrom->skip) {

      // insert a->score into array
      int j;
      for (j = 0; j < k; j++)
        if (a->score > arr[j])
          break;
      for (int m = k; m > j; m--)
        arr[m] = arr[m - 1];
      arr[j] = a->score;
      k++;

    }
  }

  *count = *count > 10 ? 10 : *count - 1; // update count of alns to keep
  *score = arr[*count - 1];               // save new min. score
}

/* int processSingle()
 * Process a set of singleton alignments, weighted to
 *   1/n (number of valid alignments).
 *   Return 1 if valid alignments found, else 0.
 */
int processSingle(char* qname, Aln* aln, int alnLen,
    bool extendOpt, int extend, bool avgExtOpt,
    Aln*** unpair, int* unpairIdx, int* unpairLen,
    int* unpairMem, float score, float asDiff,
    bool first, bool atacOpt, int atacLen5,
    int atacLen3, File bed, bool bedOpt, bool gzOut,
    bool ctrl, int sample, bool verbose) {

  // adjust AS tolerance for secondary alns
  if (score != NOSCORE)
    score -= asDiff;

  // determine number of valid single alignments
  //   (within score threshold and not to skipped chrom)
  uint8_t count = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (! a->paired && a->first == first
        && a->score >= score && a->chrom->save
        && ! a->chrom->skip)
      count++;
  }
  if (! count)
    return 0;

  // adjust score so that num alns is OK (1/2/3/4/5/6/8/10)
  if (count > 10 || count == 7 || count == 9)
    subsampleSingle(aln, alnLen, first, &count, &score);

  // find singletons to save
  uint8_t saved = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (! a->paired && a->first == first
        && a->score >= score && a->chrom->save
        && ! a->chrom->skip) {

      if (avgExtOpt)
        // for average-extension option, save alignment
        //   for later processing by processAvgExt()
        saveAvgExt(qname, a, count, unpair,
          unpairIdx, unpairLen, unpairMem);

      else
        // for other options, save singleton interval
        saveSingle(qname, a, count, extendOpt, extend,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);

      if (++saved == count)
        break;  // in case of AS ties
    }
  }

  // check for error saving alignments
  if (saved != count) {
    char* msg = (char*) memalloc(MAX_ALNS);
    sprintf(msg, "Saved %d alignments for read %s; should have been %d",
      saved, qname, count);
    exit(error(msg, ERRISSUE));
  }

  return 1;
}

/* void subsamplePair()
 * For sets of paired alns at an invalid count (>10, 9, 7),
 *   find a more stringent score.
 */
void subsamplePair(Aln* aln, int alnLen, uint8_t* count,
    float* score) {

  // insertion sort of aln scores
  float arr[*count];  // sorted array of aln scores
  int k = 0;          // count of valid alns analyzed
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired && a->full && a->score >= *score
        && a->chrom->save && ! a->chrom->skip) {

      // insert a->score into array
      int j;
      for (j = 0; j < k; j++)
        if (a->score > arr[j])
          break;
      for (int m = k; m > j; m--)
        arr[m] = arr[m - 1];
      arr[j] = a->score;
      k++;

    }
  }

  *count = *count > 10 ? 10 : *count - 1; // update count of alns to keep
  *score = arr[*count - 1];               // save new min. score
}

/* int processPair()
 * Process a set of paired alignments, weighted to
 *   1/n (number of valid alignments).
 *   Return 1 if valid alignments found, else 0.
 */
int processPair(char* qname, Aln* aln, int alnLen,
    double* totalLen, float score, float asDiff,
    bool atacOpt, int atacLen5, int atacLen3,
    File bed, bool bedOpt, bool gzOut, bool ctrl,
    int sample, bool verbose) {

  // adjust AS tolerance for secondary alns
  if (score != NOSCORE)
    score -= asDiff;

  // determine number of valid paired alignments
  //   (within score threshold and not to skipped chrom)
  uint8_t count = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired && a->full && a->score >= score
        && a->chrom->save && ! a->chrom->skip)
      count++;
  }
  if (! count)
    return 0;

  // adjust score so that num alns is OK (1/2/3/4/5/6/8/10)
  if (count > 10 || count == 7 || count == 9)
    subsamplePair(aln, alnLen, &count, &score);

  // find full fragments to save
  uint64_t fragLen = 0;     // local sum of fragment lengths
  uint8_t saved = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired && a->full && a->score >= score
        && a->chrom->save && ! a->chrom->skip) {

      // save full fragment
      fragLen += saveFragment(qname, a, count,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);

      if (++saved == count)
        break;  // in case of AS ties
    }
  }

  // check for error saving alignments
  if (saved != count) {
    char* msg = (char*) memalloc(MAX_ALNS);
    sprintf(msg, "Saved %d alignments for read %s; should have been %d",
      saved, qname, count);
    exit(error(msg, ERRISSUE));
  }

  *totalLen += (double) fragLen / count;
  return 1;
}

/* void processAlns()
 * Control processing of a set of alignments.
 *   Determine if the set has complete paired
 *   alignments or not, and what the best alignment
 *   scores are. Pass results to processPair()
 *   or processSingle().
 */
void processAlns(char* qname, Aln* aln, int alnLen,
    double* totalLen, int* pairedPr, int* singlePr,
    int* orphan, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, Aln*** unpair,
    int* unpairIdx, int* unpairLen, int* unpairMem,
    float asDiff, bool atacOpt, int atacLen5,
    int atacLen3, File bed, bool bedOpt, bool gzOut,
    bool ctrl, int sample, bool verbose) {

  // determine if paired alns are valid, and best score
  float scorePr = NOSCORE, scoreR1 = NOSCORE,
    scoreR2 = NOSCORE;
  bool pair = false;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired) {
      if (a->full) {
        // valid paired aln
        if (! pair || scorePr < a->score)
          scorePr = a->score; // best score so far
        pair = true;
      } else
        (*orphan)++;  // incomplete paired alignment
    } else if (singleOpt && ! pair) {
      // update best scores of singletons
      if (a->first && scoreR1 < a->score)
        scoreR1 = a->score;
      else if (! a->first && scoreR2 < a->score)
        scoreR2 = a->score;
    }
  }

//int prevPr = *pairedPr, prevSn = *singlePr;
  if (pair)
    // process paired alignments
    *pairedPr += processPair(qname, aln, alnLen,
      totalLen, scorePr, asDiff, atacOpt, atacLen5,
      atacLen3, bed, bedOpt, gzOut, ctrl, sample,
      verbose);

  else if (singleOpt)
    // process singleton alignments (separately for R1, R2)
    for (int j = 0; j < 2; j++)
      *singlePr += processSingle(qname, aln, alnLen,
        extendOpt, extend, avgExtOpt,
        unpair, unpairIdx, unpairLen, unpairMem,
        j ? scoreR2 : scoreR1, asDiff, ! j,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);

/*
//if (scoreR1 != NOSCORE && scoreR2 != NOSCORE) {
  if (pair)
    fprintf(stderr, "paired; score %f\n", scorePr);
  else
    fprintf(stderr, "single; scoreR1 %f, scoreR2 %f\n", scoreR1, scoreR2);
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    fprintf(stderr, "Aln %d: %s%s%s; %s:%d-%d; AS=%.0f\n",
      i, a->paired && a->full ? "paired" : "single",
      ! a->paired && a->first ? " R1" : "",
      ! a->paired && ! a->first ? " R2": "",
      a->chrom->name, a->pos[0], a->pos[1], a->score);
  }
  fprintf(stderr, "printed: %d paired, %d single\n",
    *pairedPr - prevPr, *singlePr - prevSn);
while(!getchar()) ;
}
// */

}

/*** Save alignment information ***/

/* void updatePairedAln()
 * Complete a properly paired alignment.
 */
void updatePairedAln(Aln* a, uint16_t flag,
    uint32_t pos, int length, float score) {
  if (flag & 0x40)
    a->pos[0] = flag & 0x10 ? pos + length : pos;
  else
    a->pos[1] = flag & 0x10 ? pos + length : pos;
  if (score == NOSCORE)
    a->score = NOSCORE;
  else if (a->score != NOSCORE)
    a->score += score;
  a->full = true;
}

/* bool savePairedAln()
 * Start a properly paired alignment. Return false
 *   if max. number has been reached.
 */
bool savePairedAln(Aln** aln, int* alnLen,
    uint16_t flag, Chrom* chrom, uint32_t pos,
    int length, uint32_t pnext, float score) {

  // check for excessive alignments
  if (*alnLen == MAX_ALNS)
    return false;

  // save aln info
  Aln* a = *aln + *alnLen;
  a->chrom = chrom;
  a->score = score;
  a->primary = (bool) (!(flag & 0x100));
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
  return true;
}

/* bool saveSingleAln()
 * Save the information for a singleton alignment.
 *   Return false if max. number has been reached.
 */
bool saveSingleAln(Aln** aln, int* alnLen,
    uint16_t flag, Chrom* chrom, uint32_t pos,
    int length, float score) {

  // check for excessive alignments
  if (*alnLen == MAX_ALNS)
    return false;

  // save aln info
  Aln* a = *aln + *alnLen;
  a->chrom = chrom;
  a->score = score;
  a->primary = (bool) (!(flag & 0x100));
  a->paired = false;
  a->strand = (bool) (!(flag & 0x10));
  a->first = (bool) (flag & 0x40);
  a->pos[0] = pos;
  a->pos[1] = pos + length;
  (*alnLen)++;
  return true;
}

/* bool parseAlign()
 * Parse a SAM/BAM alignment record. Save alignment
 *   info to Aln* array. Return true unless the max.
 *   number of alignments has been reached.
 */
bool parseAlign(Aln** aln, int* alnLen, uint16_t flag,
    Chrom* chrom, uint32_t pos, int length, uint32_t pnext,
    int* paired, int* single, int* secPair, int* secSingle,
    int* skipped, bool singleOpt, float score) {

  if ((flag & 0x3) == 0x3) {
    // paired alignment: save alignment information
    if (chrom->skip || ! chrom->save)
      (*skipped)++;
    else {
      (*paired)++;
      if (flag & 0x100)
        (*secPair)++;
    }

    // search for matching paired alignment (already analyzed)
    for (int i = 0; i < *alnLen; i++) {
      Aln* a = *aln + i;
      if ( a->paired && ! a->full && a->chrom == chrom
          && (flag & 0x40 ? a->pos[0] == pos : a->pos[1] == pos)
          && (flag & 0x100 ? ! a->primary : a->primary) ) {
        // complete paired alignment
        updatePairedAln(a, flag, pos, length, score);
        return true;
      }
    }

    // not found: start new paired alignment
    return savePairedAln(aln, alnLen, flag, chrom,
      pos, length, pnext, score);
  }

  // unpaired alignment
  if (chrom->skip || ! chrom->save)
    (*skipped)++;
  else {
    (*single)++;
    if (flag & 0x100)
      (*secSingle)++;
  }

  // save alignment info
  if (singleOpt)
    return saveSingleAln(aln, alnLen, flag, chrom,
      pos, length, score);

  return true;
}

/*** Save SAM/BAM header info ***/

/* void saveXBed()
 * Save BED regions to be excluded for this chrom.
 */
void saveXBed(Chrom* c, int xBedLen, Bed* xBed,
    bool verbose) {
  for (int i = 0; i < xBedLen; i++) {
    Bed* b = xBed + i;
    if (!strcmp(c->name, b->name)) {

      // check if interval is located off end of chrom
      if (b->pos[0] >= c->len) {
        if (verbose) {
          fprintf(stderr, "Warning! BED interval (%s, %d - %d) ignored\n",
            b->name, b->pos[0], b->pos[1]);
          fprintf(stderr, "  - located off end of reference (length %d)\n",
            c->len);
        }
        continue;
      }

      // insert interval into array, sorted by start pos
      int j;
      for (j = 0; j < c->bedLen; j++)
        if (b->pos[0] <= c->bedSt[j])
          break;
      c->bedLen++;
      c->bedSt = (uint32_t*) memrealloc(c->bedSt,
        c->bedLen * sizeof(uint32_t));
      c->bedEnd = (uint32_t*) memrealloc(c->bedEnd,
        c->bedLen * sizeof(uint32_t));
      for (int k = c->bedLen - 1; k > j; k--) {
        c->bedSt[k] = c->bedSt[k-1];
        c->bedEnd[k] = c->bedEnd[k-1];
      }
      c->bedSt[j] = b->pos[0];
      c->bedEnd[j] = b->pos[1];
    }
  }

  // merge overlapping intervals
  int i = 0;
  while (i < c->bedLen) {

    // check for interval past end of chrom
    if (c->bedEnd[i] > c->len) {
      if (verbose) {
        fprintf(stderr, "Warning! BED interval (%s, %d - %d) extends ",
          c->name, c->bedSt[i], c->bedEnd[i]);
        fprintf(stderr, "past end of ref.\n  - edited to (%s, %d - %d)\n",
          c->name, c->bedSt[i], c->len);
      }
      c->bedEnd[i] = c->len;
    }

    // check for overlap with previous
    if (i && c->bedSt[i] <= c->bedEnd[i-1]) {
      if (c->bedEnd[i] > c->bedEnd[i-1])
        c->bedEnd[i-1] = c->bedEnd[i];
      // shift coordinates back one
      for (int j = i; j < c->bedLen - 1; j++) {
        c->bedSt[j] = c->bedSt[j + 1];
        c->bedEnd[j] = c->bedEnd[j + 1];
      }
      c->bedLen--;
    } else
      i++;
  }
}

/* int saveChrom()
 * If chromosome (reference sequence) has not been
 *   saved yet, save it to the array. Return the index.
 */
int saveChrom(char* name, uint32_t len, int* chromLen,
    Chrom** chrom, int xcount, char** xchrList,
    int xBedLen, Bed* xBed, bool ctrl, bool verbose) {

  // determine if chrom has been saved already
  for (int i = 0; i < *chromLen; i++) {
    Chrom* c = *chrom + i;
    if (!strcmp(c->name, name)) {
      if (c->len != len)
        exit(error(c->name, ERRCHRLEN));
      if (! ctrl)
        c->save = true;
      return i;
    }
  }

  // determine if chrom should be skipped
  bool skip = false;
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
  c->save = ! ctrl; // do not save if ref in ctrl sample only
  c->diff = NULL;
  c->treat = (Pileup*) memalloc(sizeof(Pileup));
  c->treat->end = NULL;
  c->treat->cov = NULL;
  c->treatLen = 0;
  c->treatMem = 0;
  c->ctrl = (Pileup*) memalloc(sizeof(Pileup));
  c->ctrl->end = NULL;
  c->ctrl->cov = NULL;
  c->ctrlLen = 0;
  c->ctrlMem = 0;
  c->pval = NULL;
  c->pvalLen = NULL;
  c->sample = 0;
  c->qval = NULL;

  // determine if there are regions to be skipped
  c->bedSt = NULL;
  c->bedEnd = NULL;
  c->bedLen = 0;
  if (! skip)
    saveXBed(c, xBedLen, xBed, verbose);

  (*chromLen)++;
  return *chromLen - 1;
}

/* void loadChrom()
 * Save chromosome length info from a SAM header line.
 */
void loadChrom(char* line, int* chromLen, Chrom** chrom,
    int xcount, char** xchrList, int xBedLen, Bed* xBed,
    bool ctrl, bool verbose) {
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

  // remove trailing '\n'
  int i;
  for (i = 0; name[i] != '\n' && name[i] != '\0'; i++) ;
  name[i] = '\0';
  for (i = 0; len[i] != '\n' && len[i] != '\0'; i++) ;
  len[i] = '\0';

  // save chrom info to array (*chrom)
  saveChrom(name, (uint32_t) getInt(len), chromLen, chrom,
    xcount, xchrList, xBedLen, xBed, ctrl, verbose);
}

/* void checkHeader()
 * Check SAM header line for useful information:
 *   sort order or chromosome lengths.
 */
void checkHeader(char* line, int* chromLen, Chrom** chrom,
    int xcount, char** xchrList, int xBedLen, Bed* xBed,
    bool ctrl, bool verbose) {

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
    int i;
    for (i = 0; order[i] != '\n' && order[i] != '\0'; i++) ;
    order[i] = '\0';

    // sort order not unknown or coordinate
    if (order == NULL || ! strcmp(order, "unknown")
        || ! strcmp(order, "coordinate"))
      exit(error("", ERRSORT));

  } else if (! strcmp(tag, "@SQ"))
    // load chrom lengths from header line
    loadChrom(line, chromLen, chrom, xcount, xchrList,
      xBedLen, xBed, ctrl, verbose);

}

/*** SAM parsing ***/

/* bool loadFields()
 * Load alignment info from a SAM record.
 *   Return false on failure.
 */
bool loadFields(uint16_t* flag, char** rname, uint32_t* pos,
    uint8_t* mapq, char** cigar, char** rnext, uint32_t* pnext,
    int32_t* tlen, char** seq, char** qual, char** extra) {
  *extra = NULL;  // reset 'extra' fields
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
    }
    if (++i > 11) {
      *extra = strtok(NULL, "\n");
      break;
    }
    field = strtok(NULL, TAB);
  }
  return i > 11;
}

/* float getScore()
 * Search SAM optional fields for an alignment score.
 *   Return NOSCORE if not found.
 */
float getScore(char* extra) {
  if (extra == NULL)
    return NOSCORE;
  char* end;
  char* field = strtok_r(extra, TAB, &end);
  while (field != NULL) {
    char* tag = strtok(field, COL);
    if (!strcmp(tag, SCORE)) {
      for (int i = 0; i < 2; i++) {
        tag = strtok(NULL, COL);
        if (tag == NULL)
          return NOSCORE;
        if (i)
          return getFloat(tag);
      }
    }
    field = strtok_r(NULL, TAB, &end);
  }
  return NOSCORE;
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
    char* readName, double* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ, int minMapQ,
    int xcount, char** xchrList, int xBedLen, Bed* xBed,
    int* secPair, int* secSingle, int* orphan,
    int* chromLen, Chrom** chrom, bool singleOpt,
    bool extendOpt, int extend, bool avgExtOpt,
    Aln*** unpair, int* unpairMem, float asDiff,
    bool atacOpt, int atacLen5, int atacLen3, File bed,
    bool bedOpt, bool gzOut, bool ctrl, int sample,
    bool verbose) {

  // SAM fields to save
  char* qname, *rname, *cigar, *rnext, *seq, *qual, *extra;
  uint16_t flag;
  uint32_t pos, pnext;
  int32_t tlen;
  uint8_t mapq;

  int alnLen = 0;     // number of alignments for this read
  int unpairIdx = 0;  // \ indexes into singleton array(s)
  int unpairLen = 0;  // /   (with avgExtOpt)
  bool pastHeader = false;  // to check for misplaced header lines
  int count = 0;
  while (getLine(line, MAX_SIZE, in, gz) != NULL) {

    if (line[0] == '@') {
      if (pastHeader)
        exit(error(line, ERRHEAD));
      checkHeader(line, chromLen, chrom, xcount, xchrList,
        xBedLen, xBed, ctrl, verbose);
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
          extendOpt, extend, avgExtOpt, unpair,
          &unpairIdx, &unpairLen, unpairMem, asDiff,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);
      alnLen = 0;
      strncpy(readName, qname, MAX_ALNS);
    }

    // save alignment information
    int length = calcDist(qname, seq, cigar); // distance to 3' end
    float score = getScore(extra);
    if (! parseAlign(aln, &alnLen, flag, ref, pos, length,
        pnext, paired, single, secPair, secSingle, skipped,
        singleOpt, score) && verbose)
      fprintf(stderr, "Warning! Read %s has more than %d alignments\n",
        qname, MAX_ALNS);
    // NOTE: the following SAM fields are ignored:
    //   rnext, tlen, qual
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt, unpair,
      &unpairIdx, &unpairLen, unpairMem, asDiff,
      atacOpt, atacLen5, atacLen3, bed, bedOpt,
      gzOut, ctrl, sample, verbose);

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    processAvgExt(*unpair, unpairIdx, unpairLen,
      *totalLen, *pairedPr, bed, bedOpt, gzOut,
      ctrl, sample, verbose);

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

/* int arrayLen()
 * Calculate length of BAM auxiliary field of
 *   array type ('B').
 */
int arrayLen(char* extra) {
  int size = 0;
  switch (extra[0]) {
    case 'c': size = sizeof(int8_t); break;
    case 'C': size = sizeof(uint8_t); break;
    case 's': size = sizeof(int16_t); break;
    case 'S': size = sizeof(uint16_t); break;
    case 'i': size = sizeof(int32_t); break;
    case 'I': size = sizeof(uint32_t); break;
    case 'f': size = sizeof(float); break;
    default: ;
      char msg[4] = "' '";
      msg[1] = extra[0];
      exit(error(msg, ERRTYPE));
  }
  extra++;
  int32_t count = loadInt32(&extra);
  return 1 + sizeof(int32_t) + size * count;
}

/* int64_t loadInt()
 * Load an arbitrarily-sized int in little-endian format
 *   from a char* block. Return an int64_t that should
 *   be cast by the caller.
 */
int64_t loadInt(char* block, size_t size) {
  int64_t ans = 0;
  for (int i = 0; i < size; i++)
    ans = ans | ((block[i] & 0xFF) << (i*8));
  return ans;
}

/* float getBAMscore()
 * Search BAM auxiliary fields for an alignment score.
 *   Return NOSCORE if not found.
 */
float getBAMscore(char* extra, int len) {
  if (extra == NULL)
    return NOSCORE;

  // check auxiliary field
  char tag[3], val;
  tag[2] = '\0';
  int i = 0;
  while (i < len - 4) {

    // load tag and value
    for (int j = 0; j < 3; j++) {
      if (j < 2)
        tag[j] = extra[i];
      else
        val = extra[i];
      i++;
    }

    if (! strcmp(tag, SCORE)) {

      // return alignment score (cast to float)
      extra += i;
      switch (val) {
        case 'c':
          return (float) (int8_t) loadInt(extra, sizeof(int8_t));
        case 'C':
          return (float) (uint8_t) loadInt(extra, sizeof(uint8_t));
        case 's':
          return (float) (int16_t) loadInt(extra, sizeof(int16_t));
        case 'S':
          return (float) (uint16_t) loadInt(extra, sizeof(uint16_t));
        case 'i':
          return (float) (int32_t) loadInt(extra, sizeof(int32_t));
        case 'I':
          return (float) (uint32_t) loadInt(extra, sizeof(uint32_t));
        default: ;
          char msg[4] = "' '";
          msg[1] = val;
          exit(error(msg, ERRTYPE));
      }

    } else {

      // skip to next auxiliary field
      switch (val) {
        case 'A': i++; break;
        case 'c': i += sizeof(int8_t); break;
        case 'C': i += sizeof(uint8_t); break;
        case 's': i += sizeof(int16_t); break;
        case 'S': i += sizeof(uint16_t); break;
        case 'i': i += sizeof(int32_t); break;
        case 'I': i += sizeof(uint32_t); break;
        case 'f': i += sizeof(float); break;
        case 'Z': for (; extra[i] != '\0'; i++) ; i++; break;
        case 'H': for (; extra[i] != '\0'; i += 2) ; i++; break;
        case 'B': i += arrayLen(extra + i); break;
        default: ;
          char msg[4] = "' '";
          msg[1] = val;
          exit(error(msg, ERRTYPE));
      }
    }

    // check if field has gone past end of block
    if (i > len)
      exit(error("", ERRAUX));
  }

  return NOSCORE;
}

/* int parseBAM()
 * Parse the alignments in a BAM file.
 */
int parseBAM(gzFile in, char* line, Aln** aln,
    char* readName, int chromLen, Chrom* chrom, int n_ref,
    int idx[], double* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ, int minMapQ,
    int* secPair, int* secSingle, int* orphan,
    bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt, Aln*** unpair, int* unpairMem,
    float asDiff, bool atacOpt, int atacLen5,
    int atacLen3, File bed, bool bedOpt, bool gzOut,
    bool ctrl, int sample, bool verbose) {

  // BAM fields to save
  int32_t refID, pos, l_seq, next_refID, next_pos, tlen;
  uint16_t n_cigar_op, flag;
  uint8_t mapq;
  uint32_t* cigar;
  uint8_t* seq;
  char* read_name, *qual, *extra;

  int alnLen = 0;     // number of alignments for this read
  int unpairIdx = 0;  // \ indexes into singleton array(s)
  int unpairLen = 0;  // /   (with avgExtOpt)
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
          extendOpt, extend, avgExtOpt, unpair,
          &unpairIdx, &unpairLen, unpairMem, asDiff,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);
      alnLen = 0;
      strncpy(readName, read_name, MAX_ALNS);
    }

    // save alignment information
    int length = calcDistBAM(l_seq, n_cigar_op, cigar); // distance to 3' end
    float score = getBAMscore(extra, block_size - (int) (extra - line));
    if (! parseAlign(aln, &alnLen, flag, chrom + idx[refID],
        pos, length, next_pos, paired, single, secPair,
        secSingle, skipped, singleOpt, score) && verbose)
      fprintf(stderr, "Warning! Read %s has more than %d alignments\n",
        read_name, MAX_ALNS);
    // NOTE: the following BAM fields are ignored:
    //   next_refID, tlen, seq, qual
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt, unpair,
      &unpairIdx, &unpairLen, unpairMem, asDiff,
      atacOpt, atacLen5, atacLen3, bed, bedOpt,
      gzOut, ctrl, sample, verbose);

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    processAvgExt(*unpair, unpairIdx, unpairLen,
      *totalLen, *pairedPr, bed, bedOpt, gzOut,
      ctrl, sample, verbose);

  return count;
}

/* int readBAM()
 * Parse the header from a BAM file, then
 *   call parseBAM().
 */
int readBAM(gzFile in, char* line, Aln** aln,
    char* readName, double* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ, int minMapQ,
    int xcount, char** xchrList, int xBedLen, Bed* xBed,
    int* secPair, int* secSingle, int* orphan,
    int* chromLen, Chrom** chrom, bool singleOpt,
    bool extendOpt, int extend, bool avgExtOpt,
    Aln*** unpair, int* unpairMem, float asDiff,
    bool atacOpt, int atacLen5, int atacLen3, File bed,
    bool bedOpt, bool gzOut, bool ctrl, int sample,
    bool verbose) {

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
    idx[i] = saveChrom(line, (uint32_t) readInt32(in, true),
      chromLen, chrom, xcount, xchrList, xBedLen, xBed,
      ctrl, verbose);
  }

  return parseBAM(in, line, aln, readName,
    *chromLen, *chrom, n_ref, idx, totalLen, unmapped,
    paired, single, pairedPr, singlePr, supp, skipped,
    lowMapQ, minMapQ, secPair, secSingle, orphan,
    singleOpt, extendOpt, extend, avgExtOpt,
    unpair, unpairMem, asDiff, atacOpt, atacLen5,
    atacLen3, bed, bedOpt, gzOut, ctrl, sample, verbose);
}

/*** File I/O ***/

/* void openWrite()
 * Open a file for writing (stdout if file is '-'),
 *   adjusting filenames/extensions as needed.
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
  bool stdinBool = ! strcmp(inFile, "-");
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

/* int loadBED()
 * Load genomic regions to ignore from BED file.
 *   Return number saved.
 */
int loadBED(char* xFile, char* line, Bed** xBed) {
  // open BED file
  File in;
  bool gz = openRead(xFile, &in);

  // load BED records
  int count = 0;
  while (getLine(line, MAX_SIZE, in, gz) != NULL) {

    // parse BED record
    char* name = strtok(line, TAB);
    if (name == NULL)
      exit(error(line, ERRBED));
    int pos[2];
    for (int i = 0; i < 2; i++) {
      char* val = strtok(NULL, i ? TABN : TAB);
      if (val == NULL)
        exit(error(line, ERRBED));
      pos[i] = getInt(val);
    }
    if (pos[1] <= pos[0] || pos[0] < 0 || pos[1] < 0) {
      char* msg = (char*) memalloc(MAX_ALNS);
      sprintf(msg, "%s, %d - %d", name, pos[0], pos[1]);
      exit(error(msg, ERRBED));
    }

    // save info to xBed array
    *xBed = (Bed*) memrealloc(*xBed, (count + 1) * sizeof(Bed));
    Bed* b = *xBed + count;
    b->name = memalloc(1 + strlen(name));
    strcpy(b->name, name);
    b->pos[0] = pos[0];
    b->pos[1] = pos[1];
    count++;
  }

  // close file
  if ( (gz && gzclose(in.gzf) != Z_OK)
      || (! gz && fclose(in.f)) )
    exit(error(xFile, ERRCLOSE));

  return count;
}

/*** Main ***/

/* void logCounts()
 * Log alignment counts to stderr.
 */
void logCounts(int count, int unmapped, int supp,
    int skipped, Chrom* chrom, int chromLen, int minMapQ,
    int lowMapQ, int paired, int secPair, int orphan,
    int single, int secSingle, int singlePr, int pairedPr,
    double totalLen, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, bool bam, bool atacOpt,
    int atacLen) {
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
      if (c->skip || ! c->save) {
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
  if (pairedPr && ! atacOpt)
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
  if (atacOpt) {
    fprintf(stderr, "    ATAC-seq cut sites: %10d\n",
      2 * pairedPr + singlePr);
    fprintf(stderr, "      (expanded to length %dbp)\n", atacLen);
  }
}

/* void runProgram()
 * Controls the opening/closing of files,
 *   and analysis by readSAM() or readBAM().
 *   Passes results to findPeaks().
 */
void runProgram(char* inFile, char* ctrlFile, char* outFile,
    char* logFile, char* pileFile, char* bedFile,
    bool gzOut, bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt, int minMapQ, int xcount,
    char** xchrList, char* xFile, float pqvalue,
    bool qvalOpt, int minLen, int maxGap, float asDiff,
    bool atacOpt, int atacLen5, int atacLen3, bool verbose,
    int threads) {

  // open optional output files
  File bed, pile;
  if (bedFile != NULL)
    openWrite(bedFile, &bed, gzOut);
  if (pileFile != NULL)
    openWrite(pileFile, &pile, gzOut);

  // initialize variables
  char* line = (char*) memalloc(MAX_SIZE);
  int chromLen = 0;       // number of reference sequences
  Chrom* chrom = NULL;    // array of reference sequences
  Aln* aln = (Aln*) memalloc(MAX_ALNS * sizeof(Aln)); // array of saved alns
  int unpairMem = 0;      // number of unpaired alns (for avg-ext option)
  Aln** unpair = NULL;    // array of unpaired alns (for avg-ext option)
  char* readName = memalloc(MAX_ALNS + 1);  // name of read being analyzed
  readName[0] = readName[MAX_ALNS] = '\0';
  int sample = 0;         // number of sample pairs analyzed

  // save genomic regions to ignore
  Bed* xBed = NULL;
  int xBedLen = 0;
  if (xFile != NULL)
    xBedLen = loadBED(xFile, line, &xBed);

  // loop through input files (treatment and control)
  char* end1, *end2;
  char* treatName = strtok_r(inFile, COM, &end1);
  char* ctrlName = ctrlFile == NULL ? NULL
    : strtok_r(ctrlFile, COM, &end2);
  while (treatName) {

    // reset 'save' bools of each Chrom
    for (int j = 0; j < chromLen; j++)
      (chrom + j)->save = false;

    // initialize variables
    double fragLen = 0.0;   // total weighted length of all treatment fragments
    uint64_t genomeLen = 0; // reset genome length (could be dynamic)

    // process matching treatment/control files
    for (int i = 0; i < 2; i++) {

      // get treat/ctrl filename
      char* filename = treatName;
      if (i) {
        filename = ctrlName;
        if (ctrlName != NULL && !strcmp(ctrlName, "null"))
          filename = NULL;
        if (filename == NULL) {
          if (verbose)
            fprintf(stderr, "- control file #%d not provided -\n",
              sample + 1);
          savePileupNoCtrl(chrom, chromLen, fragLen,
            &genomeLen);
          break;
        }
      }

      // open input file
      File in;
      bool gz = openRead(filename, &in);
      bool bam = checkBAM(in, gz);
      if (verbose)
        fprintf(stderr, "Processing %s file #%d: %s\n",
          i ? "control" : "treatment", sample + 1,
          filename);

      // reset 'diff' array for each Chrom
      for (int j = 0; j < chromLen; j++) {
        Chrom* c = chrom + j;
        if (c->diff != NULL)
          for (int k = 0; k < 1 + c->len; k++) {
            c->diff->frac[k] = 0;
            c->diff->cov[k] = 0;
          }
      }

      // load and process alignments
      int unmapped = 0, paired = 0, single = 0, orphan = 0,
        pairedPr = 0, singlePr = 0, supp = 0, skipped = 0,
        lowMapQ = 0, secPair = 0, secSingle = 0;  // counting variables
      double totalLen = 0.0;  // total weighted length of paired fragments
      int count;
      if (bam)
        count = readBAM(in.gzf, line, &aln, readName,
          &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ,
          minMapQ, xcount, xchrList, xBedLen, xBed,
          &secPair, &secSingle, &orphan, &chromLen, &chrom,
          singleOpt, extendOpt, extend, avgExtOpt, &unpair,
          &unpairMem, asDiff, atacOpt, atacLen5, atacLen3,
          bed, bedFile != NULL, gzOut, i, sample, verbose);
      else
        count = readSAM(in, gz, line, &aln, readName,
          &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ,
          minMapQ, xcount, xchrList, xBedLen, xBed,
          &secPair, &secSingle, &orphan, &chromLen, &chrom,
          singleOpt, extendOpt, extend, avgExtOpt, &unpair,
          &unpairMem, asDiff, atacOpt, atacLen5, atacLen3,
          bed, bedFile != NULL, gzOut, i, sample, verbose);

      // log counts
      if (verbose)
        logCounts(count, unmapped, supp, skipped, chrom,
          chromLen, minMapQ, lowMapQ, paired, secPair,
          orphan, single, secSingle, singlePr, pairedPr,
          totalLen, singleOpt, extendOpt, extend,
          avgExtOpt, bam, atacOpt, atacLen5 + atacLen3);

      // save pileup values
      if (i)
        savePileupCtrl(chrom, chromLen, fragLen,
          &genomeLen, verbose);
      else
        fragLen = savePileupTreat(chrom, chromLen);

      // close input files
      if ( (gz && gzclose(in.gzf) != Z_OK)
          || (! gz && fclose(in.f)) )
        exit(error(filename, ERRCLOSE));
    }

    // calculate p-values
    savePval(chrom, chromLen, sample, pile,
      pileFile != NULL, treatName, ctrlName, gzOut);

    treatName = strtok_r(NULL, COM, &end1);
    ctrlName = ctrlFile == NULL ? NULL
      : strtok_r(NULL, COM, &end2);
    sample++;
  }

  // open output files
  File out, log;
  openWrite(outFile, &out, gzOut);
  if (logFile != NULL)
    openWrite(logFile, &log, gzOut);

  // find peaks
  findPeaks(out, log, logFile != NULL, gzOut, chrom,
    chromLen, &sample, pqvalue, qvalOpt, minLen,
    maxGap, verbose);

  // free memory
  if (xcount) {
    for (int i = 0; i < xcount; i++)
      free(xchrList[i]);
    free(xchrList);
  }
  if (xBedLen) {
    for (int i = 0; i < xBedLen; i++)
      free(xBed[i].name);
    free(xBed);
  }
  if (avgExtOpt) {
    for (int i = 0; i < unpairMem; i++)
      free(unpair[i]);
    free(unpair);
  }
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip) {
      if (chr->bedLen) {
for (int j = 0; j < chr->bedLen; j++)
  fprintf(stderr, "%s, %d - %d\n", chr->name, chr->bedSt[j], chr->bedEnd[j]);
        free(chr->bedSt);
        free(chr->bedEnd);
      }
      if (qvalOpt) {
        if (chr->qval) {
          free(chr->qval->end);
          free(chr->qval->cov);
          free(chr->qval);
        }
      }
      for (int j = 0; j < sample; j++) {
        if (chr->pval[j]) {
          free(chr->pval[j]->end);
          free(chr->pval[j]->cov);
          free(chr->pval[j]);
        }
      }
      free(chr->pval);
      free(chr->pvalLen);
      free(chr->ctrl->end);
      free(chr->ctrl->cov);
      free(chr->treat->end);
      free(chr->treat->cov);
      if (chr->diff) {
        free(chr->diff->frac);
        free(chr->diff->cov);
        free(chr->diff);
      }
    }
    free(chr->ctrl);
    free(chr->treat);
    free(chr->name);
  }
  free(chrom);
  free(aln);
  free(readName);
  free(line);

  // close files
  if ( ( gzOut && gzclose(out.gzf) != Z_OK )
      || ( ! gzOut && fclose(out.f) ) )
    exit(error(outFile, ERRCLOSE));
  if (logFile != NULL && ( ( ! gzOut && fclose(log.f) )
      || ( gzOut && gzclose(log.gzf) != Z_OK ) ) )
    exit(error(logFile, ERRCLOSE));
  if (pileFile != NULL && ( ( ! gzOut && fclose(pile.f) )
      || ( gzOut && gzclose(pile.gzf) != Z_OK ) ) )
    exit(error(pileFile, ERRCLOSE));
  if (bedFile != NULL && ( ( ! gzOut && fclose(bed.f) )
      || ( gzOut && gzclose(bed.gzf) != Z_OK ) ) )
    exit(error(bedFile, ERRCLOSE));
}

/* int saveXChrom()
 * Save list of chromosomes (ref names) to ignore.
 *   Return count.
 */
int saveXChrom(char* xchrom, char*** xchrList) {
  int i = 0;
  char* chrom = strtok(xchrom, COM);
  while (chrom != NULL) {
    *xchrList = (char**) memrealloc(*xchrList,
      (i + 1) * sizeof(char*));
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
  char* outFile = NULL, *inFile = NULL, *ctrlFile = NULL,
    *logFile = NULL, *pileFile = NULL, *bedFile = NULL,
    *xFile = NULL;
  char* xchrom = NULL;
  int extend = 0, minMapQ = 0, minLen = DEFMINLEN,
    maxGap = DEFMAXGAP, atacLen5 = DEFATAC, atacLen3 = 0,
    threads = DEFTHR;
  float asDiff = 0.0f, pqvalue = DEFQVAL;
  bool singleOpt = false, extendOpt = false,
    avgExtOpt = false, atacOpt = false, gzOut = false,
    qvalOpt = true;
  bool verbose = false;

  // parse argv
  int c;
  while ( (c = getopt_long(argc, argv, OPTIONS, long_options, NULL)) != -1 )
    switch (c) {
      case INFILE: inFile = optarg; break;
      case CTRLFILE: ctrlFile = optarg; break;
      case OUTFILE: outFile = optarg; break;
      case LOGFILE: logFile = optarg; break;
      case PILEFILE: pileFile = optarg; break;
      case BEDFILE: bedFile = optarg; break;
      case GZOPT: gzOut = true; break;
      case SINGLEOPT: singleOpt = true; break;
      case EXTENDOPT: extend = getInt(optarg); extendOpt = true; break;
      case AVGEXTOPT: avgExtOpt = true; break;
      case ATACOPT: atacOpt = true; break;
      case ATACLEN: atacLen5 = getInt(optarg); break;
      case XCHROM: xchrom = optarg; break;
      case XFILE: xFile = optarg; break;
      case MINMAPQ: minMapQ = getInt(optarg); break;
      case ASDIFF: asDiff = getFloat(optarg); break;
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
  if (atacOpt) {
    avgExtOpt = extendOpt = false;  // no singleton extensions in ATAC-seq mode
    if (atacLen5 <= 0)
      exit(error("", ERRATAC));
    // split atacLen into atacLen5 and atacLen3
    atacLen3 = (int) (atacLen5 / 2.0f + 0.5f);  // round up for 3' end
    atacLen5 /= 2;
  }
  if (asDiff < 0.0f)
    exit(error("", ERRASDIFF));
  if (threads < 1)
    exit(error("", ERRTHREAD));

  // save list of chromosomes to ignore
  int xcount = 0;
  char** xchrList = NULL;
  if (xchrom != NULL)
    xcount = saveXChrom(xchrom, &xchrList);

  // adjust significance level to -log scale
  if (pqvalue <= 0.0f || pqvalue > 1.0f)
    exit(error("", ERRPQVAL));
  pqvalue = - log10f(pqvalue);

  // send arguments to runProgram()
  runProgram(inFile, ctrlFile, outFile, logFile, pileFile,
    bedFile, gzOut, singleOpt, extendOpt, extend,
    avgExtOpt, minMapQ, xcount, xchrList, xFile, pqvalue,
    qvalOpt, minLen, maxGap, asDiff, atacOpt, atacLen5,
    atacLen3, verbose, threads);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getArgs(argc, argv);
  return 0;
}
