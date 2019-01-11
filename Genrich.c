/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.5
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
#include "Genrich.h"

/* void printVersion()
 * Print version and copyright.
 */
void printVersion(void) {
  fprintf(stderr, "Genrich, version %s\n", VERSION);
  fprintf(stderr, "Copyright (C) 2018 John M. Gaspar (jsh58@wildcats.unh.edu)\n");
  exit(EXIT_FAILURE);
}

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./Genrich  -%c <file>  -%c <file>", INFILE, OUTFILE);
  fprintf(stderr, "  [optional arguments]\n");
  fprintf(stderr, "Required arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input SAM/BAM file(s) for experimental sample(s)\n", INFILE);
  fprintf(stderr, "  -%c  <file>       Output peak file (in ENCODE narrowPeak format)\n", OUTFILE);
  fprintf(stderr, "Optional I/O arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input SAM/BAM file(s) for control sample(s)\n", CTRLFILE);
  fprintf(stderr, "  -%c  <file>       Output bedgraph-ish file for p/q values\n", LOGFILE);
  fprintf(stderr, "  -%c  <file>       Output bedgraph-ish file for pileups and p-values\n", PILEFILE);
  fprintf(stderr, "  -%c  <file>       Output BED file for reads/fragments/intervals\n", BEDFILE);
  fprintf(stderr, "  -%c  <file>       Output file for PCR duplicates (only with -%c)\n", DUPSFILE, DUPSOPT);
  fprintf(stderr, "Filtering options:\n");
  fprintf(stderr, "  -%c               Remove PCR duplicates\n", DUPSOPT);
  fprintf(stderr, "  -%c  <arg>        Comma-separated list of chromosomes to exclude\n", XCHROM);
  fprintf(stderr, "  -%c  <file>       Input BED file(s) of genomic regions to exclude\n", XFILE);
  fprintf(stderr, "  -%c  <int>        Minimum MAPQ to keep an alignment (def. 0)\n", MINMAPQ);
  fprintf(stderr, "  -%c  <float>      Keep sec alns with AS >= bestAS - <float> (def. 0)\n", ASDIFF);
  fprintf(stderr, "  -%c               Keep unpaired alignments (def. false)\n", UNPAIROPT);
  fprintf(stderr, "  -%c  <int>        Keep unpaired alns, lengths changed to <int>\n", EXTENDOPT);
  fprintf(stderr, "  -%c               Keep unpaired alns, lengths changed to paired avg\n", AVGEXTOPT);
  fprintf(stderr, "Options for ATAC-seq:\n");
  fprintf(stderr, "  -%c               Use ATAC-seq mode (def. false)\n", ATACOPT);
  fprintf(stderr, "  -%c  <int>        Expand cut sites to <int> bp (def. %d)\n", ATACLEN, DEFATAC);
  fprintf(stderr, "Options for peak-calling:\n");
  fprintf(stderr, "  -%c  <float>      Maximum q-value (FDR-adjusted p-value; def. %.2f)\n", QVALUE, DEFQVAL);
  fprintf(stderr, "  -%c  <float>      Maximum p-value (overrides -%c if set)\n", PVALUE, QVALUE);
  fprintf(stderr, "  -%c  <float>      Minimum AUC for a peak (def. %.1f)\n", MINAUC, DEFAUC);
  fprintf(stderr, "  -%c  <int>        Minimum length of a peak (def. %d)\n", MINLEN, DEFMINLEN);
  fprintf(stderr, "  -%c  <int>        Maximum distance between signif. sites (def. %d)\n", MAXGAP, DEFMAXGAP);
  fprintf(stderr, "Other options:\n");
  fprintf(stderr, "  -%c               Skip peak-calling\n", NOPEAKS);
  fprintf(stderr, "  -%c               Call peaks directly from a log file (-%c)\n", PEAKSONLY, LOGFILE);
  fprintf(stderr, "  -%c               Option to gzip-compress output(s)\n", GZOPT);
  fprintf(stderr, "  -%c               Option to print status updates/counts to stderr\n", VERBOSE);
  exit(EXIT_FAILURE);
}

/*** Utilites ***/

/* int error()
 * Prints an error message.
 */
int error(const char* msg, enum errCode err) {
  fprintf(stderr, "Error! %s%s\n", msg, errMsg[err]);
  return EXIT_FAILURE;
}

/* void* memalloc()
 * Allocates a heap block.
 */
void* memalloc(size_t size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* void* memrealloc()
 * Changes the size of a heap block.
 */
void* memrealloc(void* ptr, size_t size) {
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
    uint64_t genomeLen, float* pVal, uint64_t* pEnd,
    int64_t pLen, bool verbose) {

  // sort pileup by p-values
  quickSort(pVal, pEnd, 0, pLen - 1);

  // calculate q-values for each p-value: -log(q) = -log(p*N/k)
  uint64_t k = 1;  // 1 + number of bases with higher -log(p)
  float logN = -log10f(genomeLen);
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
      if (chr->pval[n]->cov[j] == SKIP)
        chr->qval->cov[j] = SKIP; // skipped region
      else
        chr->qval->cov[j] = lookup(pVal, 0, pLen,
          qVal, chr->pval[n]->cov[j]);
  }

  // check if all q-values are 1
  if (verbose && qVal[pLen-1] == 0.0f)
    fprintf(stderr, "Warning! All q-values are 1\n");

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
      if (p->cov[m] != SKIP)
        *pLen += recordPval(table, p->cov[m],
          p->end[m] - start);
      start = p->end[m];
    }
  }

  return table;
}

/* float* collectPval()
 * Collect arrays of p-values and genome lengths from
 *   hashtable (to be used in q-value calculations).
 */
float* collectPval(Hash** table, uint64_t** pEnd,
    int64_t pLen, uint64_t* checkLen) {
  float* pVal = (float*) memalloc(pLen * sizeof(float));
  int64_t idx = 0;
  for (int i = 0; i < HASH_SIZE; i++)
    for (Hash* h = table[i]; h != NULL; h = h->next) {
      pVal[idx] = h->val;
      (*pEnd)[idx] = h->len;
      *checkLen += h->len;
      idx++;
    }
  if (idx != pLen)
    exit(error(errMsg[ERRPVAL], ERRISSUE));
  return pVal;
}

/* void computeQval()
 * Control q-value calculations.
 */
void computeQval(Chrom* chrom, int chromLen,
    uint64_t genomeLen, int n, bool verbose) {

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
  uint64_t checkLen = 0;  // should match genomeLen
  float* pVal = collectPval(table, &pEnd, pLen, &checkLen);

  // check that collected p-value lengths match genomeLen
  if (checkLen != genomeLen) {
    char msg[MAX_ALNS];
    sprintf(msg, "Genome length (%ld) does not match p-value length (%ld)",
      genomeLen, checkLen);
    exit(error(msg, ERRISSUE));
  }

  // convert p-values to q-values
  saveQval(chrom, chromLen, n, genomeLen, pVal, pEnd,
    pLen, verbose);

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
 *   Argument 'n' is an integer in [1, 199].
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
 *   Argument 'alph' is an integer in [2, 200].
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
 *   'df' must be an even integer in [4, 400].
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
    if (pval[j] != NULL && pval[j]->cov[idx[j]] != SKIP) {
      sum += pval[j]->cov[idx[j]];
      df += 2;
    }
  if (df == 0)
    return SKIP;
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

/*** Log pileups and stats ***/

/* void printLogHeader()
 * Print header of logfile.
 */
void printLogHeader(File log, bool gzOut, int n,
    bool qvalOpt, bool sigOpt) {
  if (n) {
    // multiple samples: logfile has multiple p-values, no pileups
    if (gzOut) {
      gzprintf(log.gzf, "chr\tstart\tend");
      for (int i = 0; i < n; i++)
        gzprintf(log.gzf, "\t-log(p)_%d", i);
      gzprintf(log.gzf, "\t-log(p)_comb");
      if (qvalOpt)
        gzprintf(log.gzf, "\t-log(q)");
      if (sigOpt)
        gzprintf(log.gzf, "\tsignif");
      gzprintf(log.gzf, "\n");
    } else {
      fprintf(log.f, "chr\tstart\tend");
      for (int i = 0; i < n; i++)
        fprintf(log.f, "\t-log(p)_%d", i);
      fprintf(log.f, "\t-log(p)_comb");
      if (qvalOpt)
        fprintf(log.f, "\t-log(q)");
      if (sigOpt)
        fprintf(log.f, "\tsignif");
      fprintf(log.f, "\n");
    }
  } else {
    // single sample: logfile has pileups and p-/q-values
    if (gzOut) {
      gzprintf(log.gzf, "chr\tstart\tend\texperimental\tcontrol\t-log(p)");
      if (qvalOpt)
        gzprintf(log.gzf, "\t-log(q)");
      if (sigOpt)
        gzprintf(log.gzf, "\tsignif");
      gzprintf(log.gzf, "\n");
    } else {
      fprintf(log.f, "chr\tstart\tend\texperimental\tcontrol\t-log(p)");
      if (qvalOpt)
        fprintf(log.f, "\t-log(q)");
      if (sigOpt)
        fprintf(log.f, "\tsignif");
      fprintf(log.f, "\n");
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
    uint32_t idx[], float pval, bool qvalOpt,
    float qval, bool sig) {
  if (gzOut) {
    gzprintf(out.gzf, "%s\t%d\t%d", name, start, end);
    for (int i = 0; i < n; i++)
      if (p[i] == NULL || p[i]->cov[idx[i]] == SKIP)
        gzprintf(out.gzf, "\t%s", NA);
      else
        gzprintf(out.gzf, "\t%f", p[i]->cov[idx[i]]);
    if (pval == SKIP) {
      gzprintf(out.gzf, "\t%s", NA);
      if (qvalOpt)
        gzprintf(out.gzf, "\t%s", NA);
    } else {
      gzprintf(out.gzf, "\t%f", pval);
      if (qvalOpt)
        gzprintf(out.gzf, "\t%f", qval);
    }
    gzprintf(out.gzf, "%s\n", sig ? "\t*" : "");
  } else {
    fprintf(out.f, "%s\t%d\t%d", name, start, end);
    for (int i = 0; i < n; i++)
      if (p[i] == NULL || p[i]->cov[idx[i]] == SKIP)
        fprintf(out.f, "\t%s", NA);
      else
        fprintf(out.f, "\t%f", p[i]->cov[idx[i]]);
    if (pval == SKIP) {
      fprintf(out.f, "\t%s", NA);
      if (qvalOpt)
        fprintf(out.f, "\t%s", NA);
    } else {
      fprintf(out.f, "\t%f", pval);
      if (qvalOpt)
        fprintf(out.f, "\t%f", qval);
    }
    fprintf(out.f, "%s\n", sig ? "\t*" : "");
  }
}

/* void printInterval()
 * Print bedgraph(ish) interval for a single replicate.
 *   Values: pileups (experimental and control), -log(p),
 *   -log(q), and significance ('*') for each.
 */
void printInterval(File out, bool gzOut, char* name,
    uint32_t start, uint32_t end, float exptVal,
    float ctrlVal, float pval, bool qvalOpt, float qval,
    bool sig) {
  if (gzOut) {
    if (ctrlVal == SKIP) {
      gzprintf(out.gzf, "%s\t%d\t%d\t%f\t%f\t%s",
        name, start, end, exptVal, 0.0f, NA);
      if (qvalOpt)
        gzprintf(out.gzf, "\t%s", NA);
      gzprintf(out.gzf, "\n");
    } else {
      gzprintf(out.gzf, "%s\t%d\t%d\t%f\t%f\t%f",
        name, start, end, exptVal, ctrlVal, pval);
      if (qvalOpt)
        gzprintf(out.gzf, "\t%f", qval);
      gzprintf(out.gzf, "%s\n", sig ? "\t*" : "");
    }
  } else {
    if (ctrlVal == SKIP) {
      fprintf(out.f, "%s\t%d\t%d\t%f\t%f\t%s",
        name, start, end, exptVal, 0.0f, NA);
      if (qvalOpt)
        fprintf(out.f, "\t%s", NA);
      fprintf(out.f, "\n");
    } else {
      fprintf(out.f, "%s\t%d\t%d\t%f\t%f\t%f",
        name, start, end, exptVal, ctrlVal, pval);
      if (qvalOpt)
        fprintf(out.f, "\t%f", qval);
      fprintf(out.f, "%s\n", sig ? "\t*" : "");
    }
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
      chr->expt->cov[j], chr->ctrl->cov[k],
      chr->pval[n]->cov[m], qvalOpt,
      qvalOpt ? chr->qval->cov[m] : SKIP, sig);
  } else {
    // multiple replicates
    printIntervalN(log, gzOut, chr->name,
      start, chr->pval[n]->end[m], chr->pval, n, idx,
      chr->pval[n]->cov[m], qvalOpt,
      qvalOpt ? chr->qval->cov[m] : SKIP, sig);
    // update indexes into pval arrays
    for (int r = 0; r < n; r++)
      if (chr->pval[r] != NULL
          && chr->pval[r]->end[idx[r]] == chr->pval[n]->end[m])
        idx[r]++;
  }
}

/* void logIntervals()
 * Instead of calling peaks, just print log of pileups,
 *   and p- and q-values for each interval.
 */
void logIntervals(File log, bool gzOut, Chrom* chrom,
    int chromLen, int n, bool qvalOpt) {
  // print header
  printLogHeader(log, gzOut, n, qvalOpt, false);

  // loop through chroms
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || (qvalOpt && chr->qval == NULL)
        || (! qvalOpt && chr->pval[n] == NULL) )
      continue;

    // create indexes into arrays (expt/ctrl pileup [if n == 0]
    //   and p-value arrays [if n > 0])
    uint32_t j = 0, k = 0;  // indexes into chr->expt, chr->ctrl
    uint32_t idx[n];        // indexes into each pval array
    for (int r = 0; r < n; r++)
      idx[r] = 0;

    // loop through intervals (defined by chr->pval[n])
    uint32_t start = 0;    // start of interval
    for (uint32_t m = 0; m < chr->pvalLen[n]; m++) {

      // print stats for interval
      printLog(log, gzOut, chr, start, n, m, j, k,
        idx, qvalOpt, false);

      // update chr->expt and chr->ctrl indexes
      if (! n) {
        if (chr->ctrl->end[k] < chr->expt->end[j])
          k++;
        else {
          if (chr->ctrl->end[k] == chr->expt->end[j])
            k++;
          j++;
        }
      }

      start = chr->pval[n]->end[m];
    }
  }
}

/*** Call peaks ***/

/* void printPeak()
 * Print peaks in ENCODE narrowPeak format.
 */
void printPeak(File out, bool gzOut, char* name,
    int64_t start, int64_t end, int count, float signal,
    float pval, float qval, uint32_t pos) {
  if (gzOut) {
    gzprintf(out.gzf, "%s\t%ld\t%ld\tpeak_%d\t%d\t.\t%f\t%f",
      name, start, end, count,
      MIN((unsigned int) (1000.0f * signal / (end - start)
        + 0.5f), 1000),
      signal, pval);
    if (qval == SKIP)
      gzprintf(out.gzf, "\t-1\t%d\n", pos);
    else
      gzprintf(out.gzf, "\t%f\t%d\n", qval, pos);
  } else {
    fprintf(out.f, "%s\t%ld\t%ld\tpeak_%d\t%d\t.\t%f\t%f",
      name, start, end, count,
      MIN((unsigned int) (1000.0f * signal / (end - start)
        + 0.5f), 1000),
      signal, pval);
    if (qval == SKIP)
      fprintf(out.f, "\t-1\t%d\n", pos);
    else
      fprintf(out.f, "\t%f\t%d\n", qval, pos);
  }
}

/* void checkPeak()
 * Check if potential peak is valid (coordinates,
 *   minAUC, and minLen parameters). If valid, print
 *   results via printPeak().
 */
void checkPeak(File out, bool gzOut, char* name,
    int64_t start, int64_t end, int* count, float auc,
    float pval, float qval, uint32_t pos, float minAUC,
    int minLen, uint64_t* peakBP) {
  if (start != -1 && auc >= minAUC
      && end - start >= minLen) {
    printPeak(out, gzOut, name, start, end, *count,
      auc, pval, qval, pos);
    (*peakBP) += end - start;
    (*count)++;
  }
}

/* void resetVars()
 * Reset peak variables to null values.
 */
void resetVars(int64_t* peakStart, float* summitVal,
    uint32_t* summitLen, float* auc) {
  *peakStart = -1;
  *summitVal = -1.0f;
  *summitLen = 0;
  *auc = 0.0f;
}

/* void updatePeak()
 * Update peak variables for current interval.
 */
void updatePeak(int64_t* peakStart, int64_t* peakEnd,
    uint32_t start, uint32_t end, float* auc, float pqval,
    float minPQval, float pval, float qval,
    float* summitVal, float* summitPval, float* summitQval,
    uint32_t* summitPos, uint32_t* summitLen) {
  // update peak AUC, coordinates
  uint32_t len = end - start;
  *auc += len * (pqval - minPQval); // sum AUC
  if (*peakStart == -1)
    *peakStart = start; // start new potential peak
  *peakEnd = end;       // end of potential peak (for now)

  // check if interval is summit for this peak
  if (pqval > *summitVal) {
    *summitVal = pqval;
    *summitPval = pval;
    *summitQval = qval;
    *summitPos = (end + start)/2 - *peakStart; // midpoint of interval
    *summitLen = len;
  } else if (pqval == *summitVal) {
    // update summitPos only if interval is longer
    if (len > *summitLen) {
      *summitPos = (end + start)/2 - *peakStart; // midpoint of interval
      *summitLen = len;
      // assume summitPval, summitQval remain the same
    }
  }
}

/* int callPeaks()
 * Call peaks, using minAUC, maxGap, and minLen parameters.
 *   Produce output on the fly. Log pileups, p- and
 *   q-values for each interval. Return number of peaks.
 */
int callPeaks(File out, File log, bool logOpt, bool gzOut,
    Chrom* chrom, int chromLen, int n, float minPQval,
    bool qvalOpt, float minAUC, int minLen, int maxGap,
    uint64_t* peakBP) {

  if (logOpt)
    printLogHeader(log, gzOut, n, qvalOpt, true);

  // loop through chroms
  int count = 0;      // count of peaks
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || (qvalOpt && chr->qval == NULL)
        || (! qvalOpt && chr->pval[n] == NULL) )
      continue;

    // create indexes into arrays for logging purposes
    //   (expt/ctrl pileup [if n == 0] and p-value arrays [if n > 0])
    uint32_t j = 0, k = 0;  // indexes into chr->expt, chr->ctrl
    uint32_t idx[n];        // indexes into each pval array
    for (int r = 0; r < n; r++)
      idx[r] = 0;

    // initialize peak variables
    float auc = 0.0f;                     // area under the curve (signif.)
    int64_t peakStart = -1, peakEnd = -1; // ends of potential peak
    float summitVal = -1.0f;              // summit p/q value
    uint32_t summitPos = 0;               // distance from peakStart to summit
    uint32_t summitLen = 0;               // length of summit interval
    float summitPval = -1.0f, summitQval = -1.0f; // summit p- and q-values

    // loop through intervals (defined by chr->pval[n])
    uint32_t start = 0;    // start of interval
    for (uint32_t m = 0; m < chr->pvalLen[n]; m++) {

      bool sig = false;
      float pqval = qvalOpt ? chr->qval->cov[m]
        : chr->pval[n]->cov[m];
      if (pqval > minPQval) {

        // interval reaches significance
        sig = true;
        updatePeak(&peakStart, &peakEnd, start,
          chr->pval[n]->end[m], &auc, pqval,
          minPQval, chr->pval[n]->cov[m],
          qvalOpt ? chr->qval->cov[m] : SKIP,
          &summitVal, &summitPval, &summitQval,
          &summitPos, &summitLen);

      } else {

        // interval does not reach significance --
        //   check if interval is to be skipped
        //   OR distance is beyond maxGap from peak
        if (pqval == SKIP
            || chr->pval[n]->end[m] - peakEnd > maxGap) {
          // check if previous peak is valid
          checkPeak(out, gzOut, chr->name, peakStart,
            peakEnd, &count, auc, summitPval, summitQval,
            summitPos, minAUC, minLen, peakBP);

          // reset peak variables
          resetVars(&peakStart, &summitVal, &summitLen, &auc);
        }
      }

      // print stats for interval
      if (logOpt)
        printLog(log, gzOut, chr, start, n, m, j, k,
          idx, qvalOpt, sig);

      // update chr->expt and chr->ctrl indexes
      if (! n) {
        if (chr->ctrl->end[k] < chr->expt->end[j])
          k++;
        else {
          if (chr->ctrl->end[k] == chr->expt->end[j])
            k++;
          j++;
        }
      }

      start = chr->pval[n]->end[m];
    }

    // determine if last peak is valid
    checkPeak(out, gzOut, chr->name, peakStart, peakEnd,
      &count, auc, summitPval, summitQval, summitPos,
      minAUC, minLen, peakBP);
  }

  return count;
}

/* void findPeaks()
 * Control process of finding peaks:
 *   calculating p- and q-values, calling peaks,
 *   and printing output.
 */
void findPeaks(File out, File log, bool logOpt, bool gzOut,
    Chrom* chrom, int chromLen, int* sample,
    float minPQval, bool qvalOpt, int minLen, int maxGap,
    float minAUC, bool peaksOpt, bool verbose) {

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
      for (int j = 0; j < chr->bedLen; j += 2)
        genomeLen -= chr->bed[j+1] - chr->bed[j];
    }
  }

  if (verbose) {
    if (peaksOpt) {
      fprintf(stderr, "Peak-calling parameters:\n");
      fprintf(stderr, "  Genome length: %ldbp\n", genomeLen);
      fprintf(stderr, "  Significance threshold: -log(%c) > %.3f\n",
        qvalOpt ? 'q' : 'p', minPQval);
      fprintf(stderr, "  Min. AUC: %.3f\n", minAUC);
      if (minLen)
        fprintf(stderr, "  Min. peak length: %dbp\n", minLen);
      fprintf(stderr, "  Max. gap between sites: %dbp\n", maxGap);
    } else {
      fprintf(stderr, "- peak-calling skipped -\n");
      fprintf(stderr, "  Genome length: %ldbp\n", genomeLen);
    }
  }

  // compute q-values
  if (qvalOpt)
    computeQval(chrom, chromLen, genomeLen, *sample - 1,
      verbose);

  // call peaks
  if (peaksOpt) {
    uint64_t peakBP = 0;
    int count = callPeaks(out, log, logOpt, gzOut, chrom,
      chromLen, *sample - 1, minPQval, qvalOpt, minAUC,
      minLen, maxGap, &peakBP);
    if (verbose)
      fprintf(stderr, "Peaks identified: %d (%ldbp)\n",
        count, peakBP);
  } else if (logOpt)
    // not calling peaks -- just log pileups and stats
    logIntervals(log, gzOut, chrom, chromLen, *sample - 1,
      qvalOpt);
}

/*** Call peaks directly from a log file ***/

/* void saveXBed()
 * Save BED regions to be excluded for a chrom.
 */
void saveXBed(char* name, uint32_t len, int* bedLen,
    uint32_t** bed, int xBedLen, Bed* xBed,
    bool verbose) {
  // check each xBed interval for match to chrom (name)
  for (int i = 0; i < xBedLen; i++) {
    Bed* b = xBed + i;
    if (!strcmp(name, b->name)) {

      // check if interval is located off end of chrom
      if (b->pos[0] >= len) {
        if (verbose) {
          fprintf(stderr, "Warning! BED interval (%s, %d - %d) ignored\n",
            b->name, b->pos[0], b->pos[1]);
          fprintf(stderr, "  - located off end of reference %s (length %d)\n",
            name, len);
        }
        continue;
      }

      // insert interval into array, sorted by start pos
      int j;
      for (j = 0; j < *bedLen; j += 2)
        if (b->pos[0] <= (*bed)[j])
          break;
      *bedLen += 2;
      *bed = (uint32_t*) memrealloc(*bed,
        *bedLen * sizeof(uint32_t));
      // shift intervals forward
      for (int k = *bedLen - 1; k > j + 1; k--)
        (*bed)[k] = (*bed)[k-2];
      (*bed)[j] = b->pos[0];
      (*bed)[j+1] = b->pos[1];
    }
  }

  // merge overlapping intervals
  int i = 0;
  while (i < *bedLen) {

    // check for interval past end of chrom
    if ((*bed)[i+1] > len) {
      if (verbose) {
        fprintf(stderr, "Warning! BED interval (%s, %d - %d) extends ",
          name, (*bed)[i], (*bed)[i+1]);
        fprintf(stderr, "past end of ref.\n  - edited to (%s, %d - %d)\n",
          name, (*bed)[i], len);
      }
      (*bed)[i+1] = len;
    }

    // check for overlap with previous
    if (i && (*bed)[i] <= (*bed)[i-1]) {
      if ((*bed)[i+1] > (*bed)[i-1])
        (*bed)[i-1] = (*bed)[i+1];
      // shift coordinates backward
      for (int j = i; j < *bedLen - 2; j++)
        (*bed)[j] = (*bed)[j + 2];
      *bedLen -= 2;
    } else
      i += 2;

  }
}

/* bool checkChrom()
 * Return true if given char* matches a chromosome
 *   to be skipped (-e).
 */
bool checkChrom(char* chr, int xcount, char** xchrList) {
  for (int i = 0; i < xcount; i++)
    if (!strcmp(xchrList[i], chr))
      return true;
  return false;
}

/* int getIdx()
 * Find fields of header line of bedgraph-ish log file
 *   that match "-log([pq])".
 *   If multiple matches, return the last one.
 */
int getIdx(File in, bool gz, char* line, bool qvalOpt,
    int* idxQ) {
  if (getLine(line, MAX_SIZE, in, gz) == NULL)
    exit(error("<header>", ERRLOGIDX));
  char matchP[8] = "-log(p)";
  char matchQ[8] = "-log(q)";
  int idxP = -1;
  int i = 0;
  char* field = strtok(line, TABN);
  while (field != NULL) {
    if (!strncmp(field, matchP, 7))
      idxP = i;
    else if (!strncmp(field, matchQ, 7))
      *idxQ = i;
    field = strtok(NULL, TABN);
    i++;
  }
  if (idxP == -1)
    exit(error(matchP, ERRLOGIDX));
  if (qvalOpt && *idxQ == -1)
    exit(error(matchQ, ERRLOGIDX));
  return idxP;
}

/* void loadBDG()
 * Load fields (chr, start, end, and [pq]-value) from
 *   a bedgraph-ish record (-f log file).
 */
void loadBDG(char* line, char** chr, uint32_t* start,
    uint32_t* end, char** pStat, int idxP, char** qStat,
    int idxQ, bool qvalOpt) {
  int idx = qvalOpt ? idxQ : idxP;
  char* field = strtok(line, TABN);
  for (int i = 0; i <= idx; i++) {
    if (field == NULL)
      exit(error("", ERRLOG));
    switch (i) {
      case CHR: *chr = field; break;
      case START: *start = getInt(field); break;
      case END: *end = getInt(field); break;
      default:
        if (i == idxP)
          *pStat = field;
        else if (i == idxQ)
          *qStat = field;
    }
    field = strtok(NULL, TABN);
  }
}

/* void callPeaksLog()
 * Call peaks directly from a bedgraph-ish log file.
 */
void callPeaksLog(File in, bool gz, File out, bool gzOut,
    char* line, int xcount, char** xchrList, int xBedLen,
    Bed* xBed, float minPQval, bool qvalOpt, int minLen,
    int maxGap, float minAUC, bool verbose) {
  // determine index of fields ([pq]-value) to analyze
  int idxQ = -1;
  int idxP = getIdx(in, gz, line, qvalOpt, &idxQ);

  // initialize variables for genomic regions to be excluded
  uint32_t* bed = NULL;         // array of BED coordinates
  int bedLen = 0;               // length of bed array
  int bedIdx = 0;               // index into bed array
  uint32_t bedPos = UINT32_MAX; // next coordinate
  bool save = true;             // status of current interval
  bool warn = false;            // warn of new skipped regions?

  // initialize counting variables
  int count = 0;
  uint64_t genomeLen = 0;
  uint64_t peakBP = 0;

  // initialize peak variables
  float auc = 0.0f;                     // area under the curve (signif.)
  int64_t peakStart = -1, peakEnd = -1; // ends of potential peak
  float summitVal = -1.0f;              // summit p/q value
  uint32_t summitPos = 0;               // distance from peakStart to summit
  uint32_t summitLen = 0;               // length of summit interval
  float summitPval = -1.0f, summitQval = -1.0f; // summit p- and q-values

  // parse bedgraph records
  char prev[MAX_ALNS];              // previous chrom name
  prev[0] = '\0';
  bool skip = false;                // skip current chrom?
  char* chr, *stat, *pStat, *qStat; // chrom and stats/NA of current interval
  uint32_t start, end;              // coordinates of current interval
  float pqval;                      // [pq]-value of current interval
  while (getLine(line, MAX_SIZE, in, gz) != NULL) {

    // load fields from bedgraph record
    loadBDG(line, &chr, &start, &end, &pStat, idxP,
      &qStat, idxQ, qvalOpt);

    // new chromosome
    if (strcmp(prev, chr)) {
      // check if previous chrom's last peak is valid
      checkPeak(out, gzOut, prev, peakStart, peakEnd,
        &count, auc, summitPval, summitQval, summitPos,
        minAUC, minLen, &peakBP);

      // reset peak variables
      resetVars(&peakStart, &summitVal, &summitLen, &auc);

      // check if new chrom should be skipped
      skip = checkChrom(chr, xcount, xchrList);
      if (verbose && skip) {
        fprintf(stderr, "Warning! Skipping chromosome %s --\n  ", chr);
        fprintf(stderr, "Reads aligning to it were used in the background");
        fprintf(stderr, " pileup calculation,\n  and its length was included");
        fprintf(stderr, " in the genome length %scalculation\n",
          qvalOpt ? "(and q-value) " : "");
      }

      // check for regions to be skipped
      bedLen = 0;
      if (! skip) {
        // load regions from xBed
        // (note: cannot check validity of coordinates)
        saveXBed(chr, UINT32_MAX, &bedLen, &bed, xBedLen,
          xBed, verbose);

        // load next coordinate
        bedIdx = 0;
        bedPos = bedIdx < bedLen ? bed[bedIdx]
          : UINT32_MAX;
        save = true;
      }

      strcpy(prev, chr);  // update current chrom
    }
    if (skip)
      continue; // current chrom to be skipped

    // check stat for 'NA' (skipped region)
    stat = qvalOpt ? qStat : pStat;
    if (!strcmp(stat, NA)) {
      // check if previous peak is valid
      checkPeak(out, gzOut, chr, peakStart, peakEnd,
        &count, auc, summitPval, summitQval, summitPos,
        minAUC, minLen, &peakBP);

      // reset peak variables
      resetVars(&peakStart, &summitVal, &summitLen, &auc);
      continue;
    }
    pqval = getFloat(stat);

    // check if current interval starts at a position
    //   to be skipped (new -E)
    if (bedPos == start) {
      if (save) {
        // check if previous peak is valid
        checkPeak(out, gzOut, chr, peakStart, peakEnd,
          &count, auc, summitPval, summitQval, summitPos,
          minAUC, minLen, &peakBP);
        // reset peak variables
        resetVars(&peakStart, &summitVal, &summitLen, &auc);
      }

      // load next coordinate
      save = ! save;
      bedIdx++;
      bedPos = bedIdx < bedLen ? bed[bedIdx] : UINT32_MAX;
    }

    // check if current interval contains the next position
    //   to be skipped (new -E)
    uint32_t subStart = start;  // start of current subinterval
    while (bedPos > start && bedPos < end) {

      if (save) {
        // interval *not* to be skipped:
        //   update peak if necessary, then print if valid
        if (pqval > minPQval) {
          // interval reaches significance
          updatePeak(&peakStart, &peakEnd, subStart,
            bedPos, &auc, pqval, minPQval,
            qvalOpt ? getFloat(pStat) : pqval,
            qvalOpt ? pqval : SKIP,
            &summitVal, &summitPval, &summitQval,
            &summitPos, &summitLen);
        }
        // check if peak is valid
        checkPeak(out, gzOut, chr, peakStart, peakEnd,
          &count, auc, summitPval, summitQval, summitPos,
          minAUC, minLen, &peakBP);
        // reset peak variables
        resetVars(&peakStart, &summitVal, &summitLen, &auc);
        genomeLen += bedPos - subStart; // update genome length
      } else
        warn = true;

      // load next coordinate
      subStart = bedPos;
      save = ! save;
      bedIdx++;
      bedPos = bedIdx < bedLen ? bed[bedIdx] : UINT32_MAX;
    }
    if (! save) {
      warn = true;
      continue;
    }
    start = subStart; // reset start coordinate

    // check [pq]-value for significance
    genomeLen += end - start; // update genome length
    if (pqval > minPQval) {

      // interval reaches significance
      updatePeak(&peakStart, &peakEnd, start, end,
        &auc, pqval, minPQval,
        qvalOpt ? getFloat(pStat) : pqval,
        qvalOpt ? pqval : SKIP,
        &summitVal, &summitPval, &summitQval,
        &summitPos, &summitLen);

    } else {

      // interval does not reach significance --
      //   check if distance is beyond maxGap from peak
      if (end - peakEnd > maxGap) {
        // check if previous peak is valid
        checkPeak(out, gzOut, chr, peakStart, peakEnd,
          &count, auc, summitPval, summitQval, summitPos,
          minAUC, minLen, &peakBP);

        // reset peak variables
        resetVars(&peakStart, &summitVal, &summitLen, &auc);
      }
    }

  }

  // check if last peak is valid
  checkPeak(out, gzOut, chr, peakStart, peakEnd,
    &count, auc, summitPval, summitQval, summitPos,
    minAUC, minLen, &peakBP);

  if (verbose) {
    if (warn) {
      fprintf(stderr, "Warning! Skipping given BED regions --\n  ");
      fprintf(stderr, "Reads aligning to them were used in the background");
      fprintf(stderr, " pileup calculation,\n  and the lengths were included");
      fprintf(stderr, " in the genome length %scalculation\n",
        qvalOpt ? "(and q-value) " : "");
    }
    fprintf(stderr, "Peak-calling parameters:\n");
    fprintf(stderr, "  Genome length: %ldbp\n", genomeLen);
    fprintf(stderr, "  Significance threshold: -log(%c) > %.3f\n",
      qvalOpt ? 'q' : 'p', minPQval);
    fprintf(stderr, "  Min. AUC: %.3f\n", minAUC);
    if (minLen)
      fprintf(stderr, "  Min. peak length: %dbp\n", minLen);
    fprintf(stderr, "  Max. gap between sites: %dbp\n", maxGap);
    fprintf(stderr, "Peaks identified: %d (%ldbp)\n",
      count, peakBP);
  }

  free(bed);
}

/*** Calculate a p-value (log-normal distribution) ***/
// adapted from R-3.5.0 source code, as noted below

/* double do_del()
 * Adapted from pnorm.c in R-3.5.0
 *   (cf. do_del() and swap_tail).
 */
double do_del(double y, double temp, bool ret) {
  double xsq = trunc(y * 16) / 16;
  double del = (y - xsq) * (y + xsq);
  if (ret)
    return log1p(-exp((-xsq * xsq - del) / 2.0) * temp);
  return (-xsq * xsq - del) / 2.0 + log(temp);
}

/* double pnorm()
 * Adapted from pnorm.c in R-3.5.0
 *   (cf. pnorm_both() with i_tail=1, log_p=TRUE).
 */
double pnorm(double x) {
  double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };
  double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };
  double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  double xden, xnum, xsq, temp;
  double y = fabs(x);
  if (y <= 0.67448975) {

    // small values of fabs(x)
    if (y > DBL_EPSILON * 0.5) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (int i = 0; i < 3; i++) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
      temp = x * (xnum + a[3]) / (xden + b[3]);
    } else
      temp = x * a[3] / b[3];
    return log(0.5 - temp);

  } else if (y <= sqrt(32.0)) {

    // slightly larger values of fabs(x)
    xnum = c[8] * y;
    xden = y;
    for (int i = 0; i < 7; i++) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);
    return do_del(y, temp, x <= 0.0);

  } else if (y < 1e170) {

    // even larger values of fabs(x)
    xsq = 1.0 / (x * x);
    xnum = p[5] * xsq;
    xden = xsq;
    for (int i = 0; i < 4; i++) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (1/sqrt(2*M_PI) - temp) / y;
    return do_del(x, temp, x <= 0.0);
  }

  // default
  return -0.0;
}

/* double plnorm()
 * Calculate a p-value for a log-normal distribution
 *   with observation 'x' and parameters 'meanlog' and
 *   'sdlog'.
 * Adapted from plnorm.c and pnorm.c in R-3.5.0,
 *   with lower_tail=FALSE and log_p=TRUE.
 * Return value is -log10(p).
 */
double plnorm(double x, double meanlog, double sdlog) {
  if (sdlog == 0.0)
    return x < meanlog ? 0.0 : FLT_MAX;
  return -pnorm((log(x) - meanlog) / sdlog) / M_LN10;
}

/* float calcPval()
 * Calculate -log10(p) using a log-normal distribution
 *   with mu=ctrlVal, sd={mu>7 ? 10*log10(mu) : 1.2*mu},
 *   and observation exptVal.
 */
float calcPval(float exptVal, float ctrlVal) {
  if (ctrlVal == SKIP)
    return SKIP; // in a skipped region
  if (ctrlVal == 0.0f)
    return exptVal == 0.0f ? 0.0f : FLT_MAX;
  if (exptVal == 0.0f)
    return 0.0f;

  // calculate meanlog and sdlog for plnorm()
  double meanlog, sdlog;
  double mu = ctrlVal;
  if (mu > 7.0) {
    double sd = 10.0 * log10(mu);
    mu *= mu;
    sd *= sd;
    meanlog = log(mu / sqrt(sd + mu));
    sdlog = sqrt(log1p(sd / mu));
  } else {
    meanlog = log(mu) - LOGSQRT;
    sdlog = SQRTLOG;
  }

  // calculate pval by plnorm()
  double pval = plnorm(exptVal, meanlog, sdlog);
  return pval > FLT_MAX ? FLT_MAX : (float) pval;
}

/*** Create and save p-values genome-wide ***/

/* uint32_t countIntervals()
 * Count the number of pileup intervals to create
 *   for a composite.
 */
uint32_t countIntervals(Chrom* chr) {
  uint32_t num = 0;
  uint32_t k = 0;
  for (uint32_t j = 0; j < chr->exptLen; j++) {
    while (k < chr->ctrlLen
        && chr->ctrl->end[k] < chr->expt->end[j]) {
      num++;
      k++;
    }
    if (chr->ctrl->end[k] == chr->expt->end[j])
      k++;
    num++;
  }
  return num;
}

/* void printPileHeader()
 * Print header of the bedgraph-ish pileup log file.
 */
void printPileHeader(File pile, char* exptName,
    char* ctrlName, bool gzOut) {
  if (gzOut) {
    gzprintf(pile.gzf, "# experimental file: %s; control file: %s\n",
      exptName, ctrlName && strcmp(ctrlName, "null") ? ctrlName : NA);
    gzprintf(pile.gzf, "chr\tstart\tend\texperimental\tcontrol\t-log(p)\n");
  } else {
    fprintf(pile.f, "# experimental file: %s; control file: %s\n",
      exptName, ctrlName && strcmp(ctrlName, "null") ? ctrlName : NA);
    fprintf(pile.f, "chr\tstart\tend\texperimental\tcontrol\t-log(p)\n");
  }
}

/* void printPile()
 * Print bedgraph-ish interval of experimental/control
 *   pileup values and p-value.
 */
void printPile(File pile, char* name, uint32_t start,
    uint32_t end, float expt, float ctrl, float pval,
    bool gzOut) {
  if (gzOut) {
    if (ctrl == SKIP)
      gzprintf(pile.gzf, "%s\t%d\t%d\t%f\t%f\t%s\n",
        name, start, end, expt, 0.0f, NA);
    else
      gzprintf(pile.gzf, "%s\t%d\t%d\t%f\t%f\t%f\n",
        name, start, end, expt, ctrl, pval);
  } else {
    if (ctrl == SKIP)
      fprintf(pile.f, "%s\t%d\t%d\t%f\t%f\t%s\n",
        name, start, end, expt, 0.0f, NA);
    else
      fprintf(pile.f, "%s\t%d\t%d\t%f\t%f\t%f\n",
        name, start, end, expt, ctrl, pval);
  }
}

/* void savePval()
 * Create and save p-values as pileups for each Chrom*.
 */
void savePval(Chrom* chrom, int chromLen, int n,
    File pile, bool pileOpt, char* exptName,
    char* ctrlName, bool gzOut) {

  // print log header
  if (pileOpt)
    printPileHeader(pile, exptName, ctrlName, gzOut);

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
      if (chr->ctrl->end[k] < chr->expt->end[j]) {
        p->end[m] = chr->ctrl->end[k];
        p->cov[m] = calcPval(chr->expt->cov[j],
          chr->ctrl->cov[k]);
        if (pileOpt)
          printPile(pile, chr->name, start, p->end[m],
            chr->expt->cov[j], chr->ctrl->cov[k],
            p->cov[m], gzOut);
        k++;
      } else {
        p->end[m] = chr->expt->end[j];
        p->cov[m] = calcPval(chr->expt->cov[j],
          chr->ctrl->cov[k]);
        if (pileOpt)
          printPile(pile, chr->name, start, p->end[m],
            chr->expt->cov[j], chr->ctrl->cov[k],
            p->cov[m], gzOut);
        if (chr->ctrl->end[k] == chr->expt->end[j])
          k++;
        j++;
      }
      start = p->end[m];
    }

  }
}

/*** Save experimental/control pileup values ***/

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
    double fragLen) {
  uint64_t genomeLen = 0;
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip && chr->save) {
      genomeLen += chr->len;
      for (int j = 0; j < chr->bedLen; j += 2)
        genomeLen -= chr->bed[j+1] - chr->bed[j];
    }
  }
  if (! genomeLen)
    exit(error("", ERRGEN));
  return fragLen / genomeLen;
}

/* void saveLambda()
 * For a ctrl pileup, save constant value of lambda,
 *   or -1 (SKIP) for BED intervals to be skipped.
 */
void saveLambda(Chrom* chr, float lambda) {
  if (chr->bedLen == 0) {
    // no BED intervals: save constant of lambda
    saveConst(chr->ctrl, &chr->ctrlLen, &chr->ctrlMem,
      chr->len, lambda);
    return;
  }

  // with BED intervals
  int num = chr->bedLen + 1;  // number of array intervals
  int idx = 0;    // index into chr->bed array
  bool save = true;
  if (chr->bed[0] == 0) {
    num--;
    idx++;
    save = false; // 1st interval skipped
  }
  if (chr->bed[chr->bedLen - 1] == chr->len)
    num--;

  // expand pileup arrays (if necessary)
  if (num > chr->ctrlMem) {
    chr->ctrl->end = (uint32_t*) memrealloc(chr->ctrl->end,
      num * sizeof(uint32_t));
    chr->ctrl->cov = (float*) memrealloc(chr->ctrl->cov,
      num * sizeof(float));
    chr->ctrlMem = num;
  }
  chr->ctrlLen = num;

  // populate chr->ctrl arrays: alternate lambda and -1
  for (int j = 0; j < num - 1; j++) {
    chr->ctrl->end[j] = chr->bed[idx];
    chr->ctrl->cov[j] = save ? lambda : SKIP;
    save = ! save;
    idx++;
  }
  chr->ctrl->end[num-1] = chr->len;
  chr->ctrl->cov[num-1] = save ? lambda : SKIP;
}

/* void savePileupNoCtrl()
 * When no control is available, save the control
 *   pileup as the background lambda value.
 */
void savePileupNoCtrl(Chrom* chrom, int chromLen,
    double fragLen, bool verbose) {
  float lambda = calcLambda(chrom, chromLen, fragLen);
  if (verbose)
    fprintf(stderr, "  Background pileup value: %f\n",
      lambda);
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save)
      continue;
    saveLambda(chr, lambda);
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

  return getVal(*cov, *frac);
}

/* float calcFactor
 * Calculate the scaling factor of experimental fragment
 *   lengths to control. Also, set ctrlLen for each
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

    // initialize BED interval variables
    int bedIdx = 0;         // index into chr->bed array
    uint32_t bedPos = bedIdx < chr->bedLen ?
      chr->bed[bedIdx] : chr->len + 1;  // position of BED interval
    bool save = true;
    if (bedPos == 0) {
      // first BED interval starts at 0
      save = false;
      bedIdx++;
      bedPos = bedIdx < chr->bedLen ?
        chr->bed[bedIdx] : chr->len + 1;
    }

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

      if (j == bedPos || (save && (d->cov[j] || d->frac[j]))) {
        if (save)
          // save frag. length weighted by val
          ctrlFrag += (j - start) * val;
        num++;
        start = j;
      }

      // update pileup value
      if (d->cov[j] || d->frac[j])
        val = updateVal(d->cov[j], d->frac[j], &cov, &frac);

      // update 'save' status (from BED intervals)
      if (j == bedPos) {
        save = ! save;
        bedIdx++;
        bedPos = bedIdx < chr->bedLen ?
          chr->bed[bedIdx] : chr->len + 1;
      }
    }

    // save final interval
    if (save)
      ctrlFrag += (j - start) * val;
    chr->ctrlLen = num; // save number of intervals
  }

  // return ratio of experimental frags to ctrl frags
  if (! ctrlFrag)
    return 1.0f;
  return fragLen / ctrlFrag;
}

/* void savePileupCtrl()
 * Save pileup values for a control sample from
 *   'diff' arrays and background lambda value.
 */
void savePileupCtrl(Chrom* chrom, int chromLen,
    double fragLen, bool verbose) {

  // calculate background lambda value
  float lambda = calcLambda(chrom, chromLen, fragLen);
  if (verbose)
    fprintf(stderr, "  Background pileup value: %f\n", lambda);

  // calculate scale factor (experimental / control)
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
      saveLambda(chr, lambda);
      continue;
    }

    // expand pileup arrays (if necessary)
    //   (chr->ctrlLen already set in calcFactor())
    if (chr->ctrlLen > chr->ctrlMem) {
      chr->ctrl->end = (uint32_t*) memrealloc(chr->ctrl->end,
        chr->ctrlLen * sizeof(uint32_t));
      chr->ctrl->cov = (float*) memrealloc(chr->ctrl->cov,
        chr->ctrlLen * sizeof(float));
      chr->ctrlMem = chr->ctrlLen;
    }

    // initialize BED interval variables
    int bedIdx = 0;         // index into chr->bed array
    uint32_t bedPos = bedIdx < chr->bedLen ?
      chr->bed[bedIdx] : chr->len + 1;  // position of BED interval
    bool save = true;
    if (bedPos == 0) {
      // first BED interval starts at 0 (should be skipped)
      save = false;
      bedIdx++;
      bedPos = bedIdx < chr->bedLen ?
        chr->bed[bedIdx] : chr->len + 1;
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

      // update pileup value
      if (d->cov[j] || d->frac[j])
        val = factor * updateVal(d->cov[j], d->frac[j],
          &cov, &frac);

      // determine if interval should be saved
      if (j == bedPos || (save && net != MAX(val, lambda))) {
        chr->ctrl->end[pos] = j;
        chr->ctrl->cov[pos] = save ? net : SKIP;
        pos++;
      }
      net = MAX(val, lambda);

      // update 'save' status (from BED intervals)
      if (j == bedPos) {
        save = ! save;
        bedIdx++;
        bedPos = bedIdx < chr->bedLen ?
          chr->bed[bedIdx] : chr->len + 1;
      }

    }

    // save final interval
    chr->ctrl->end[pos] = j;
    chr->ctrl->cov[pos] = save ? net : SKIP;

    // update array length
    if (pos >= chr->ctrlLen) {
      char msg[MAX_ALNS];
      sprintf(msg, "%s (%s)", errMsg[ERRARRC], chr->name);
      exit(error(msg, ERRISSUE));
    }
    chr->ctrlLen = pos + 1;

    // check for val error (should end at 0)
    val = updateVal(d->cov[j], d->frac[j], &cov, &frac);
    if (val) {
      char msg[MAX_ALNS];
      sprintf(msg, "Control pileup for ref %s finishes at %f (not 0.0)",
        chr->name, val);
      exit(error(msg, ERRISSUE));
    }
  }

}

/* double savePileupExpt()
 * Save pileup values for an experimental sample from
 *   'diff' arrays.
 *   Return total length of all fragments (weighted).
 */
double savePileupExpt(Chrom* chrom, int chromLen) {

  // create pileup for each chrom
  double fragLen = 0.0;  // weighted fragment length
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->skip || ! chr->save)
      continue;

    // if no read coverage, save constant pileup of 0
    if (chr->diff == NULL) {
      saveConst(chr->expt, &chr->exptLen, &chr->exptMem,
        chr->len, 0.0f);
      continue;
    }

    // determine number of pileup intervals
    int bedIdx = 0;     // index into chr->bed array
    uint32_t bedPos = bedIdx < chr->bedLen ?
      chr->bed[bedIdx] : chr->len + 1;  // position of BED interval
    bool save = true;
    if (bedPos == 0) {
      // first BED interval starts at 0 (should be skipped)
      save = false;
      bedIdx++;
      bedPos = bedIdx < chr->bedLen ?
        chr->bed[bedIdx] : chr->len + 1;
    }
    Diff* d = chr->diff;
    uint32_t num = 1;   // number of intervals
    for (uint32_t j = 1; j < chr->len; j++)
      if (j == bedPos) {
        num++;
        save = ! save;
        bedIdx++;
        bedPos = bedIdx < chr->bedLen ?
          chr->bed[bedIdx] : chr->len + 1;
      } else if (save && (d->cov[j] || d->frac[j]))
        num++;

    // expand pileup arrays (if necessary)
    if (num > chr->exptMem) {
      chr->expt->end = (uint32_t*) memrealloc(chr->expt->end,
        num * sizeof(uint32_t));
      chr->expt->cov = (float*) memrealloc(chr->expt->cov,
        num * sizeof(float));
      chr->exptMem = num;
    }
    chr->exptLen = num;

    // reset BED interval values
    bedIdx = 0; // index into chr->bed array
    bedPos = bedIdx < chr->bedLen ?
      chr->bed[bedIdx] : chr->len + 1;
    save = true;
    if (bedPos == 0) {
      save = false;
      bedIdx++;
      bedPos = bedIdx < chr->bedLen ?
        chr->bed[bedIdx] : chr->len + 1;
    }

    // initialize pileup values
    int32_t cov = 0;      // current pileup value
    uint8_t frac = 0;     // current pileup value (fraction part)
    float val = updateVal(d->cov[0], d->frac[0], &cov, &frac);

    // save pileup values along the chrom
    uint32_t start = 0;   // beginning coordinate of interval
    uint32_t pos = 0;     // position in pileup arrays
    uint32_t j;
    for (j = 1; j < chr->len; j++) {

      if (j == bedPos || (save && (d->cov[j] || d->frac[j]))) {
        // save end of interval and pileup value
        chr->expt->end[pos] = j;
        if (save) {
          chr->expt->cov[pos] = val;
          fragLen += (j - start) * val; // frag. length weighted by val
        } else
          chr->expt->cov[pos] = 0.0f;
        pos++;
        start = j;
      }

      // update pileup value
      if (d->cov[j] || d->frac[j])
        val = updateVal(d->cov[j], d->frac[j], &cov, &frac);

      // update 'save' status (from BED intervals)
      if (j == bedPos) {
        save = ! save;
        bedIdx++;
        bedPos = bedIdx < chr->bedLen ?
          chr->bed[bedIdx] : chr->len + 1;
      }

    }

    // save final interval
    chr->expt->end[pos] = j;
    if (save) {
      chr->expt->cov[pos] = val;
      fragLen += (j - start) * val;
    } else
      chr->expt->cov[pos] = 0.0f;

    // verify array length
    if (pos + 1 != chr->exptLen) {
      char msg[MAX_ALNS];
      sprintf(msg, "%s (%s)", errMsg[ERRARR], chr->name);
      exit(error(msg, ERRISSUE));
    }

    // check for val error (should end at 0)
    val = updateVal(d->cov[j], d->frac[j], &cov, &frac);
    if (val) {
      char msg[MAX_ALNS];
      sprintf(msg, "Experimental pileup for ref %s finishes at %f (not 0.0)",
        chr->name, val);
      exit(error(msg, ERRISSUE));
    }
  }

  if (fragLen == 0.0)
    exit(error("", ERREXPT));
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
      char msg[MAX_ALNS];
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
      char msg[MAX_ALNS];
      sprintf(msg, "%s (%d)", errMsg[ERRALNS], count);
      exit(error(msg, ERRISSUE));
  }
}

/*** Convert alignments to intervals ***/

/* void printBED()
 * Print a BED interval for a read/fragment.
 *   Append the aln count, 'C'ontrol/'E'xperimental,
 *   and sample number to the read name (4th column).
 */
void printBED(File bed, bool gzOut, char* chr,
    int64_t start, int64_t end, char* qname,
    uint8_t count, bool ctrl, int sample) {
  if (gzOut)
    gzprintf(bed.gzf, "%s\t%ld\t%ld\t%s_%d_%c_%d\n",
      chr, start, end, qname, count, ctrl ? 'C' : 'E',
      sample);
  else
    fprintf(bed.f, "%s\t%ld\t%ld\t%s_%d_%c_%d\n",
      chr, start, end, qname, count, ctrl ? 'C' : 'E',
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
    char msg[MAX_ALNS];
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
  if (c->diff->cov[start] == INT16_MAX) {
    if (verbose) {
      fprintf(stderr, "Warning! Read %s, alignment at (%s, %ld-%ld)",
        qname, c->name, start, end);
      fprintf(stderr, " skipped due to overflow\n");
    }
    return 0;
  }
  if (c->diff->cov[end] == INT16_MIN) {
    if (verbose) {
      fprintf(stderr, "Warning! Read %s, alignment at (%s, %ld-%ld)",
        qname, c->name, start, end);
      fprintf(stderr, " skipped due to underflow\n");
    }
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

/* int calcAvgLen()
 * Calculate the average fragment length.
 *   Return 0 if no fragments.
 */
int calcAvgLen(double totalLen, int pairedPr,
    bool verbose) {
  if (! pairedPr) {
    if (verbose) {
      fprintf(stderr, "Warning! No paired alignments to calculate avg frag ");
      fprintf(stderr, "length --\n  Printing unpaired alignments \"as is\"\n");
    }
    return 0;
  }
  return (int) (totalLen / pairedPr + 0.5);
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
  int avgLen = calcAvgLen(totalLen, pairedPr, verbose);

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
 * Save info for an unpaired alignment to an array,
 *   for "extend to average length" option. Alignments
 *   will be processed later by processAvgExt().
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

/* void saveUnpair()
 * Control processing of unpaired alignments
 *   (either keeping them as is, or extending
 *   to a given length).
 */
void saveUnpair(char* qname, Aln* a, uint8_t count,
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

/*** Save alignments for evaluation of PCR duplicates ***/

/* Read* createRead()
 * Expand read arrays if necessary and update
 *   indexes. Return pointer to next Read*.
 */
Read* createRead(Read*** read, int* readIdx,
    int* readLen, int* readMem) {
  // alloc memory if necessary
  if (*readLen == 0 && *readIdx == *readMem) {
    // check if max. number of reads exceeded
    if ((*readMem + 1) * MAX_SIZE > UINT32_MAX) {
      char msg[MAX_ALNS];
      sprintf(msg, "Exceeded max. number of reads (%u)",
        UINT32_MAX);
      exit(error(msg, ERRISSUE));
    }

    *read = (Read**) memrealloc(*read,
      (*readMem + 1) * sizeof(Read*));
    (*read)[*readMem] = (Read*) memalloc(MAX_SIZE
      * sizeof(Read));
    (*readMem)++;
  }
  Read* r = (*read)[*readIdx] + *readLen;

  // update indexes
  (*readLen)++;
  if (*readLen == MAX_SIZE) {
    *readLen = 0;
    (*readIdx)++;
  }
  return r;
}

/* void copyAlns()
 * Copy alignment info for a set of singleton
 *   alignments.
 */
void copyAlns(Aln* aln, int alnLen, float score,
    float asDiff, bool first, Aln** dest,
    uint8_t* destLen) {
  // adjust AS tolerance for secondary alns
  if (score != NOSCORE)
    score -= asDiff;

  // determine number of valid single alignments
  uint8_t count = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (! a->paired && a->first == first
        && a->score >= score)
      count++;
  }
  *dest = (Aln*) memalloc(count * sizeof(Aln));
  *destLen = count;

  // copy alignment info for valid alignments
  uint8_t j = 0;  // index into r->aln
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (! a->paired && a->first == first
        && a->score >= score) {
      Aln* b = *dest + j;
      b->paired = a->paired;
      b->first = a->first;
      b->strand = a->strand;
      b->score = a->score;
      b->chrom = a->chrom;
      b->pos[0] = a->pos[0];
      b->pos[1] = a->pos[1];
      j++;
    }
  }

}

/* void saveAlnsSingle()
 * Save a set of singleton alignments.
 */
void saveAlnsSingle(char* qname, Aln* aln, int alnLen,
    float score, float asDiff, bool first, Read* r,
    uint16_t qual) {
  // populate Read* struct
  r->name = (char*) memalloc(1 + strlen(qname));
  strcpy(r->name, qname);
  r->qual = qual;
  r->first = first;
  r->score = score;

  // copy alignments to struct
  copyAlns(aln, alnLen, score, asDiff, first, &r->aln,
    &r->alnLen);
}

/* void saveAlnsDiscord()
 * Save a set of discordant alignments.
 */
void saveAlnsDiscord(char* qname, Aln* aln, int alnLen,
    float scoreR1, float scoreR2, float asDiff,
    Read* r, uint16_t qualR1, uint16_t qualR2) {
  // save R1 alignments
  saveAlnsSingle(qname, aln, alnLen, scoreR1, asDiff,
    true, r, qualR1);
  // save R2 alignments
  copyAlns(aln, alnLen, scoreR2, asDiff, false, &r->alnR2,
    &r->alnLenR2);
  r->qual = MIN(qualR1 + qualR2, UINT16_MAX);
  r->scoreR2 = scoreR2;
}

/* void saveAlnsPair()
 * Save a set of properly paired alignments.
 */
void saveAlnsPair(char* qname, Aln* aln, int alnLen,
    float score, float asDiff, Read* r, uint16_t qualR1,
    uint16_t qualR2) {
  // populate Read* struct
  r->name = (char*) memalloc(1 + strlen(qname));
  strcpy(r->name, qname);
  r->qual = MIN(qualR1 + qualR2, UINT16_MAX);
  r->score = score;

  // adjust AS tolerance for secondary alns
  if (score != NOSCORE)
    score -= asDiff;

  // determine number of valid paired alignments
  uint8_t count = 0;
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired && a->full && a->score >= score)
      count++;
  }
  r->aln = (Aln*) memalloc(count * sizeof(Aln));
  r->alnLen = count;

  // copy alignment info for valid alignments
  uint8_t j = 0;  // index into r->aln
  for (int i = 0; i < alnLen; i++) {
    Aln* a = aln + i;
    if (a->paired && a->full && a->score >= score) {
      Aln* b = r->aln + j;
      b->paired = a->paired;
      b->full = a->full;
      b->score = a->score;
      b->chrom = a->chrom;
      // ensure positions are ordered
      if (a->pos[0] > a->pos[1]) {
        b->pos[0] = a->pos[1];
        b->pos[1] = a->pos[0];
      } else {
        b->pos[0] = a->pos[0];
        b->pos[1] = a->pos[1];
      }
      j++;
    }
  }

}

/* void saveAlns()
 * Control saving of alignments. Use createRead()
 *   to make a Read*, and pass that to saveAlnsPair(),
 *   saveAlnsDiscord(), or saveAlnsSingle().
 */
void saveAlns(char* qname, Aln* aln, int alnLen, bool pair,
    bool singleOpt, bool singleR1, bool singleR2,
    float scorePr, float scoreR1, float scoreR2,
    float asDiff, Read*** readPr, int* readIdxPr,
    int* readLenPr, int* readMemPr, Read*** readDc,
    int* readIdxDc, int* readLenDc, int* readMemDc,
    Read*** readSn, int* readIdxSn, int* readLenSn,
    int* readMemSn, uint16_t qualR1, uint16_t qualR2) {
  if (pair) {
    // properly paired alignment(s)
    Read* r = createRead(readPr, readIdxPr, readLenPr,
      readMemPr);
    saveAlnsPair(qname, aln, alnLen, scorePr, asDiff, r,
      qualR1, qualR2);
  } else if (singleOpt) {
    if (singleR1 && singleR2) {
      // both reads aligned (discordant)
      Read* r = createRead(readDc, readIdxDc, readLenDc,
        readMemDc);
      saveAlnsDiscord(qname, aln, alnLen, scoreR1,
        scoreR2, asDiff, r, qualR1, qualR2);
    } else if (singleR1) {
      // only R1 read aligned
      Read* r = createRead(readSn, readIdxSn, readLenSn,
        readMemSn);
      saveAlnsSingle(qname, aln, alnLen, scoreR1, asDiff,
        true, r, qualR1);
    } else if (singleR2) {
      // only R2 read aligned
      Read* r = createRead(readSn, readIdxSn, readLenSn,
        readMemSn);
      saveAlnsSingle(qname, aln, alnLen, scoreR2, asDiff,
        false, r, qualR2);
    }
  }
}

/*** Process a set of alignments ***/

/* void subsampleSingle()
 * For sets of unpaired alns at an invalid count (>10, 9, 7),
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
 * Process a set of unpaired alignments, weighted to
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

  // determine number of valid unpaired alignments
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

  // find unpaired alns to save
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
        // for other options, save interval
        saveUnpair(qname, a, count, extendOpt, extend,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);

      if (++saved == count)
        break;  // in case of AS ties
    }
  }

  // check for error saving alignments
  if (saved != count) {
    char msg[MAX_ALNS];
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
    char msg[MAX_ALNS];
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
 *   scores are.
 * If duplicate removal is required, save alignments via
 *   saveAlns(), else pass results to processPair()
 *   or processSingle() directly.
 */
void processAlns(char* qname, Aln* aln, int alnLen,
    double* totalLen, int* pairedPr, int* singlePr,
    int* orphan, bool singleOpt, bool extendOpt,
    int extend, bool avgExtOpt, Aln*** unpair,
    int* unpairIdx, int* unpairLen, int* unpairMem,
    float asDiff, bool atacOpt, int atacLen5,
    int atacLen3, File bed, bool bedOpt, bool gzOut,
    bool ctrl, int sample, bool dupsOpt, Read*** readPr,
    int* readIdxPr, int* readLenPr, int* readMemPr,
    Read*** readDc, int* readIdxDc, int* readLenDc,
    int* readMemDc, Read*** readSn, int* readIdxSn,
    int* readLenSn, int* readMemSn, uint16_t qualR1,
    uint16_t qualR2, bool verbose) {

  // determine if paired alns are valid, and best score
  float scorePr = NOSCORE, scoreR1 = NOSCORE,
    scoreR2 = NOSCORE;
  bool pair = false, singleR1 = false, singleR2 = false;
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
      // update best scores of unpaired alns
      if (a->first && scoreR1 < a->score) {
        scoreR1 = a->score;
        singleR1 = true;
      } else if (! a->first && scoreR2 < a->score) {
        scoreR2 = a->score;
        singleR2 = true;
      }
    }
  }

  if (dupsOpt)
    // save alignments for later evaluation of duplicates
    saveAlns(qname, aln, alnLen, pair, singleOpt, singleR1,
      singleR2, scorePr, scoreR1, scoreR2, asDiff,
      readPr, readIdxPr, readLenPr, readMemPr,
      readDc, readIdxDc, readLenDc, readMemDc,
      readSn, readIdxSn, readLenSn, readMemSn,
      qualR1, qualR2);

  else {
    // process alns directly
    if (pair)
      // process paired alignments
      *pairedPr += processPair(qname, aln, alnLen,
        totalLen, scorePr, asDiff, atacOpt, atacLen5,
        atacLen3, bed, bedOpt, gzOut, ctrl, sample,
        verbose);
    else if (singleOpt) {
      // process unpaired alignments (separately for R1, R2)
      if (singleR1)
        *singlePr += processSingle(qname, aln, alnLen,
          extendOpt, extend, avgExtOpt,
          unpair, unpairIdx, unpairLen, unpairMem,
          scoreR1, asDiff, true,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);
      if (singleR2)
        *singlePr += processSingle(qname, aln, alnLen,
          extendOpt, extend, avgExtOpt,
          unpair, unpairIdx, unpairLen, unpairMem,
          scoreR2, asDiff, false,
          atacOpt, atacLen5, atacLen3, bed, bedOpt,
          gzOut, ctrl, sample, verbose);
    }
  }
}

/*** Efficient read sorting ***/

/* uint32_t johnPartition()
 * Partition the reads of a section of the qual/order
 *   arrays based on one value into 3 bins (greater,
 *   equal, and lower).
 */
uint32_t johnPartition(uint16_t* qual, uint32_t* order,
    uint32_t low, uint32_t high, uint16_t* qual0,
    uint16_t* qual1, uint16_t* qual2, uint32_t* order0,
    uint32_t* order1, uint32_t* order2,
    uint32_t* idxHigh) {

  // separate qual values into temp arrays --
  //   qual0 for higher values, qual1 equal, qual2 lower
  uint16_t pivot = qual[high - 1];  // pivot value: last elt
  uint32_t idx0 = 0, idx1 = 0, idx2 = 0; // indexes into temp arrays
  for (uint32_t j = low; j < high; j++) {
    if (qual[j] > pivot) {
      qual0[idx0] = qual[j];
      order0[idx0] = order[j];
      idx0++;
    } else if (qual[j] == pivot) {
      qual1[idx1] = qual[j];
      order1[idx1] = order[j];
      idx1++;
    } else {
      qual2[idx2] = qual[j];
      order2[idx2] = order[j];
      idx2++;
    }
  }

  if (! idx0 && ! idx2)
    return 0; // all equal values, no need to shuffle

  // recombine temp arrays back into qual/order
  uint32_t i = 0;
  bool val0 = (bool) idx0, val1 = true;
  for (uint32_t j = low; j < high; j++) {
    if (val0) {
      qual[j] = qual0[i];
      order[j] = order0[i];
      if (++i == idx0) {
        val0 = false;
        i = 0;
      }
    } else if (val1) {
      qual[j] = qual1[i];
      order[j] = order1[i];
      if (++i == idx1) {
        val1 = false;
        i = 0;
      }
    } else {
      qual[j] = qual2[i];
      order[j] = order2[i];
      i++;
    }
  }

  // return low/high indexes
  *idxHigh = low + idx0 + idx1;
  return low + idx0;
}

/* void johnSort()
 * Variation of quicksort that is stable and optimized
 *   for repeated qual values.
 */
void johnSort(uint16_t* qual, uint32_t* order,
    uint32_t low, uint32_t high, uint16_t* qual0,
    uint16_t* qual1, uint16_t* qual2, uint32_t* order0,
    uint32_t* order1, uint32_t* order2) {
  if (low + 1 < high) {
    uint32_t idx1 = 0; // new low index for upper recursive call
    uint32_t idx = johnPartition(qual, order, low, high,
      qual0, qual1, qual2, order0, order1, order2, &idx1);
    if (idx)
      // lower recursive call
      johnSort(qual, order, low, idx, qual0, qual1,
        qual2, order0, order1, order2);
    if (idx1)
      // upper recursive call
      johnSort(qual, order, idx1, high, qual0, qual1,
        qual2, order0, order1, order2);
  }
}

/* void sortReads()
 * Determine sort order of a Read** array, based on
 *   sums of quality scores.
 * Sorting performed by johnSort(), an efficient,
 *   *stable* quicksort that is optimized for these
 *   arrays in which values are frequently repeated.
 */
void sortReads(Read** arr, uint32_t count, uint32_t* order,
    uint16_t* qual, uint32_t* order0, uint32_t* order1,
    uint32_t* order2, uint16_t* qual0, uint16_t* qual1,
    uint16_t* qual2) {

  // initialize order and qual arrays
  for (uint32_t i = 0; i < count; i++) {
    order[i] = i;
    qual[i] = (arr[i / MAX_SIZE] + i % MAX_SIZE)->qual;
  }

  // initialize johnSort()
  johnSort(qual, order, 0, count, qual0, qual1, qual2,
    order0, order1, order2);
}

/*** PCR duplicate removal ***/

/* uint32_t calcHashSize()
 * Calculate a size for a hashtable:
 *   a power of 2 >= the given count * 4/3.
 * The given count is the number of reads, which
 *   will not necessarily be the same as the number
 *   of alignments. Some reads will have multiple
 *   alignments, but they will be counterbalanced
 *   by PCR duplicates (which will be discarded).
 */
uint32_t calcHashSize(uint32_t count) {
  uint32_t size = 2;
  uint32_t val = MIN(UINT32_MAX, 4 * count / 3);
  for (int i = 1; i < 31; i++) {
    if (size >= val)
      return size;
    size *= 2;
  }
  // max size 2^31
  return size;
}

/* uint32_t jenkins_hash_aln()
 * Adapted from http://www.burtleburtle.net/bob/hash/doobs.html
 *   Modified to take an alignment as input, hashed
 *   differently depending on alignment type.
 *   Returns index into hashtable.
 */
uint32_t jenkins_hash_aln(Chrom* chrom, Chrom* chrom1,
    uint32_t pos, uint32_t pos1, bool strand, bool strand1,
    int alignType, uint32_t hashSize) {
  uint32_t hash = 0;
  unsigned char* p;

  // hash Chrom*
  int end = (alignType == DISCORD ? 2 : 1);
  for (int j = 0; j < end; j++) {
    p = (unsigned char*) (j ? chrom1 : chrom);
    for (int i = 0; i < sizeof(Chrom*); i++) {
      hash += p[i];
      hash += hash << 10;
      hash ^= hash >> 6;
    }
  }

  // hash pos
  end = (alignType == SINGLE ? 1 : 2);
  for (int j = 0; j < end; j++) {
    p = (unsigned char*) (j ? &pos1 : &pos);
    for (int i = 0; i < sizeof(uint32_t); i++) {
      hash += p[i];
      hash += hash << 10;
      hash ^= hash >> 6;
    }
  }

  // hash strand
  end = alignType;  // convenient!
  for (int j = 0; j < end; j++) {
    p = (unsigned char*) (j ? &strand1 : &strand);
    for (int i = 0; i < sizeof(bool); i++) {
      hash += p[i];
      hash += hash << 10;
      hash ^= hash >> 6;
    }
  }

  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;
  return hash % hashSize;
}

/* void addToHash()
 * Add a new node (HashAln) to hashtable, inserted
 *   at given idx (already calculated).
 */
void addToHash(Chrom* chrom, Chrom* chrom1, uint32_t pos,
    uint32_t pos1, bool strand, bool strand1,
    HashAln** table, uint32_t idx, char* name) {
  HashAln* h = (HashAln*) memalloc(sizeof(HashAln));
  h->chrom = chrom;
  h->chrom1 = chrom1;
  h->pos = pos;
  h->pos1 = pos1;
  h->strand = strand;
  h->strand1 = strand1;
  h->name = NULL;
  if (name) {
    h->name = (char*) memalloc(1 + strlen(name));
    strcpy(h->name, name);
  }
  h->next = table[idx];
  table[idx] = h;
}

/* HashAln* checkHash()
 * Check hashtable for a match to an alignment.
 *   The alignment attributes are defined by the
 *   alignment type.
 */
HashAln* checkHash(Chrom* chrom, Chrom* chrom1,
    uint32_t pos, uint32_t pos1, bool strand, bool strand1,
    int alignType, HashAln** table, uint32_t idx) {
  for (HashAln* h = table[idx]; h != NULL; h = h->next) {
    // check for match, based on alignment type
    switch (alignType) {
      case PAIRED:
        if (chrom == h->chrom && pos == h->pos
            && pos1 == h->pos1)
          return h;
        break;
      case SINGLE:
        if (chrom == h->chrom && pos == h->pos
            && strand == h->strand)
          return h;
        break;
      case DISCORD:
        if (chrom == h->chrom && chrom1 == h->chrom1
            && pos == h->pos && pos1 == h->pos1
            && strand == h->strand && strand1 == h->strand1)
          return h;
        break;
      default:
        exit(error("", ERRALNTYPE));
    }
  }
  return NULL;
}

/* void checkAndAdd()
 * Check a singleton alignment for a match to the
 *   hashtable; if there is none, add it to the table.
 */
void checkAndAdd(HashAln** tableSn, uint32_t hashSizeSn,
    Chrom* chrom, uint32_t pos, bool strand, char* name) {
  uint32_t idx = jenkins_hash_aln(chrom, NULL, pos, 0,
    strand, 0, SINGLE, hashSizeSn);
  if (! checkHash(chrom, NULL, pos, 0, strand, 0, SINGLE,
      tableSn, idx))
    addToHash(chrom, NULL, pos, 0, strand, 0, tableSn,
      idx, name);
}

/* void logDup()
 * Print log information about a read identified as a
 *   duplicate, based on alignment type.
 */
void logDup(File dups, bool gzOut, char* name,
    Chrom* chrom, Chrom* chrom1, uint32_t pos,
    uint32_t pos1, bool strand, bool strand1,
    char* match, int alignType) {
  switch (alignType) {
    case PAIRED:
      if (gzOut)
        gzprintf(dups.gzf, "%s\t%s:%d-%d\t%s\tpaired\n",
          name, chrom->name, pos, pos1, match);
      else
        fprintf(dups.f, "%s\t%s:%d-%d\t%s\tpaired\n",
          name, chrom->name, pos, pos1, match);
      break;
    case SINGLE:
      if (gzOut)
        gzprintf(dups.gzf, "%s\t%s:%d,%c\t%s\tsingle\n",
          name, chrom->name, pos, strand ? '+' : '-',
          match);
      else
        fprintf(dups.f, "%s\t%s:%d,%c\t%s\tsingle\n",
          name, chrom->name, pos, strand ? '+' : '-',
          match);
      break;
    case DISCORD:
      if (gzOut)
        gzprintf(dups.gzf, "%s\t%s:%d,%c;%s:%d,%c\t%s\tdiscordant\n",
          name, chrom->name, pos, strand ? '+' : '-',
          chrom1->name, pos1, strand1 ? '+' : '-', match);
      else
        fprintf(dups.f, "%s\t%s:%d,%c;%s:%d,%c\t%s\tdiscordant\n",
          name, chrom->name, pos, strand ? '+' : '-',
          chrom1->name, pos1, strand1 ? '+' : '-', match);
      break;
    default:
      exit(error("", ERRALNTYPE));
  }
}

/* void addHashPr()
 * Add all paired alignments for a Read* to the hashtable.
 */
void addHashPr(Read* r, HashAln** table,
    uint32_t hashSize, bool dupsVerb,
    HashAln** tableSn, uint32_t hashSizeSn) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t idx = jenkins_hash_aln(a->chrom, NULL,
      a->pos[0], a->pos[1], 0, 0, PAIRED, hashSize);
    addToHash(a->chrom, NULL, a->pos[0], a->pos[1],
      0, 0, table, idx, dupsVerb ? r->name : NULL);

    // also add both alignments as singletons to hashtable
    if (tableSn != NULL) {
      checkAndAdd(tableSn, hashSizeSn, a->chrom, a->pos[0],
        true, dupsVerb ? r->name : NULL);
      checkAndAdd(tableSn, hashSizeSn, a->chrom, a->pos[1],
        false, dupsVerb ? r->name : NULL);
    }
  }
}

/* bool checkHashPr()
 * Check a set of paired alignments for a match in the
 *   hashtable. Return true if *any* match.
 */
bool checkHashPr(Read* r, HashAln** table,
    uint32_t hashSize, File dups, bool dupsVerb,
    bool gzOut) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t idx = jenkins_hash_aln(a->chrom, NULL,
      a->pos[0], a->pos[1], 0, 0, PAIRED, hashSize);
    HashAln* h = checkHash(a->chrom, NULL, a->pos[0],
      a->pos[1], 0, 0, PAIRED, table, idx);
    if (h) {
      if (dupsVerb)
        logDup(dups, gzOut, r->name, a->chrom, NULL,
          a->pos[0], a->pos[1], 0, 0, h->name, PAIRED);
      return true;
    }
  }
  return false;
}

/* void findDupsPr()
 * Find PCR duplicates among properly paired
 *   alignment sets.
 */
void findDupsPr(Read** readPr, int readIdxPr,
    int readLenPr, int* countPr, int* dupsPr,
    int* pairedPr, double* totalLen, float asDiff,
    bool atacOpt, int atacLen5, int atacLen3, File bed,
    bool bedOpt, bool gzOut, bool ctrl, int sample,
    HashAln*** table, uint32_t* tableMem,
    HashAln** tableSn, uint32_t hashSizeSn, File dups,
    bool dupsVerb, uint32_t* order, uint32_t* order0,
    uint32_t* order1, uint32_t* order2, uint16_t* qual,
    uint16_t* qual0, uint16_t* qual1, uint16_t* qual2,
    bool verbose) {

  // initialize hashtable
  uint32_t count = readIdxPr * MAX_SIZE + readLenPr;
  uint32_t hashSize = calcHashSize(count);
  if (hashSize > *tableMem) {
    *table = (HashAln**) memrealloc(*table,
      hashSize * sizeof(HashAln*));
    *tableMem = hashSize;
  } else
    hashSize = *tableMem;
  for (uint32_t i = 0; i < hashSize; i++)
    (*table)[i] = NULL;

  // get sort order of reads by qual score sum
  sortReads(readPr, count, order, qual, order0, order1,
    order2, qual0, qual1, qual2);

  // loop through paired reads
  for (uint32_t i = 0; i < count; i++) {
    Read* r = readPr[order[i] / MAX_SIZE]
      + order[i] % MAX_SIZE;

    // check hashtable for matches
    if (checkHashPr(r, *table, hashSize, dups,
        dupsVerb, gzOut))
      (*dupsPr)++;
    else {
      // add alignments to hashtable(s)
      addHashPr(r, *table, hashSize, dupsVerb,
        tableSn, hashSizeSn);
      // process alignments too
      *pairedPr += processPair(r->name, r->aln,
        r->alnLen, totalLen, r->score, asDiff,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
    }

    // free Read
    free(r->name);
    free(r->aln);

    (*countPr)++;
  }

  // free nodes of hashtable
  for (uint32_t i = 0; i < hashSize; i++) {
    HashAln* tmp;
    HashAln* h = (*table)[i];
    while (h != NULL) {
      free(h->name);
      tmp = h->next;
      free(h);
      h = tmp;
    }
  }
}

/* void addHashDc()
 * Add all combinations of discordant alignments for a
 *   Read* to the hashtable.
 */
void addHashDc(Read* r, HashAln** table,
    uint32_t hashSize, bool dupsVerb,
    HashAln** tableSn, uint32_t hashSizeSn) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t pos = (a->strand ? a->pos[0] : a->pos[1]);
    for (int j = 0; j < r->alnLenR2; j++) {
      Aln* b = r->alnR2 + j;
      uint32_t pos1 = (b->strand ? b->pos[0] : b->pos[1]);
      uint32_t idx = jenkins_hash_aln(a->chrom, b->chrom,
        pos, pos1, a->strand, b->strand, DISCORD, hashSize);
      addToHash(a->chrom, b->chrom, pos, pos1, a->strand,
        b->strand, table, idx, dupsVerb ? r->name : NULL);

      // also add both alignments as singletons to hashtable
      if (tableSn != NULL) {
        if (! j)
          checkAndAdd(tableSn, hashSizeSn, a->chrom, pos,
            a->strand, dupsVerb ? r->name : NULL);
        if (! k)
          checkAndAdd(tableSn, hashSizeSn, b->chrom, pos1,
            b->strand, dupsVerb ? r->name : NULL);
      }
    }
  }
}

/* bool checkHashDc()
 * Check each combination of discordant alignments for a
 *   match in the hashtable. Return true if *any* match.
 */
bool checkHashDc(Read* r, HashAln** table,
    uint32_t hashSize, File dups, bool dupsVerb,
    bool gzOut) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t pos = (a->strand ? a->pos[0] : a->pos[1]);
    for (int j = 0; j < r->alnLenR2; j++) {
      Aln* b = r->alnR2 + j;
      uint32_t pos1 = (b->strand ? b->pos[0] : b->pos[1]);
      uint32_t idx = jenkins_hash_aln(a->chrom, b->chrom,
        pos, pos1, a->strand, b->strand, DISCORD, hashSize);
      HashAln* h = checkHash(a->chrom, b->chrom, pos, pos1,
        a->strand, b->strand, DISCORD, table, idx);
      if (h) {
        if (dupsVerb)
          logDup(dups, gzOut, r->name, a->chrom, b->chrom,
            pos, pos1, a->strand, b->strand, h->name,
            DISCORD);
        return true;
      }
      // check the reverse also
      idx = jenkins_hash_aln(b->chrom, a->chrom,
        pos1, pos, b->strand, a->strand, DISCORD, hashSize);
      h = checkHash(b->chrom, a->chrom, pos1, pos,
        b->strand, a->strand, DISCORD, table, idx);
      if (h) {
        if (dupsVerb)
          logDup(dups, gzOut, r->name, b->chrom, a->chrom,
            pos1, pos, b->strand, a->strand, h->name,
            DISCORD);
        return true;
      }
    }
  }
  return false;
}

/* void findDupsDc()
 * Find PCR duplicates among discordant
 *   alignment sets.
 */
void findDupsDc(Read** readDc, int readIdxDc,
    int readLenDc, int* countDc, int* dupsDc,
    int* singlePr, bool extendOpt, int extend,
    float asDiff, bool atacOpt, int atacLen5, int atacLen3,
    File bed, bool bedOpt, bool gzOut, bool ctrl,
    int sample, HashAln*** table, uint32_t* tableMem,
    HashAln** tableSn, uint32_t hashSizeSn, File dups,
    bool dupsVerb, uint32_t* order, uint32_t* order0,
    uint32_t* order1, uint32_t* order2, uint16_t* qual,
    uint16_t* qual0, uint16_t* qual1, uint16_t* qual2,
    bool verbose) {

  // initialize hashtable
  uint32_t count = readIdxDc * MAX_SIZE + readLenDc;
  uint32_t hashSize = calcHashSize(count);
  if (hashSize > *tableMem) {
    *table = (HashAln**) memrealloc(*table,
      hashSize * sizeof(HashAln*));
    *tableMem = hashSize;
  } else
    hashSize = *tableMem;
  for (uint32_t i = 0; i < hashSize; i++)
    (*table)[i] = NULL;

  // sort reads by qual score sum
  sortReads(readDc, count, order, qual, order0, order1,
    order2, qual0, qual1, qual2);

  // loop through discordant reads
  for (uint32_t i = 0; i < count; i++) {
    Read* r = readDc[order[i] / MAX_SIZE]
      + order[i] % MAX_SIZE;

    // check hashtable for matches
    if (checkHashDc(r, *table, hashSize, dups,
        dupsVerb, gzOut))
      (*dupsDc)++;
    else {
      // add alignments to hashtable(s)
      addHashDc(r, *table, hashSize, dupsVerb,
        tableSn, hashSizeSn);
      // process alignments too (as singletons)
      *singlePr += processSingle(r->name, r->aln,
        r->alnLen, extendOpt, extend, false,
        NULL, NULL, NULL, NULL,
        r->score, asDiff, true,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
      *singlePr += processSingle(r->name, r->alnR2,
        r->alnLenR2, extendOpt, extend, false,
        NULL, NULL, NULL, NULL,
        r->scoreR2, asDiff, false,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
    }

    // free Read
    free(r->name);
    free(r->aln);
    free(r->alnR2);

    (*countDc)++;
  }

  // free nodes of hashtable
  for (uint32_t i = 0; i < hashSize; i++) {
    HashAln* tmp;
    HashAln* h = (*table)[i];
    while (h != NULL) {
      free(h->name);
      tmp = h->next;
      free(h);
      h = tmp;
    }
  }
}

/* void addHashSn()
 * Add all singleton alignments for a Read* to the
 *   hashtable.
 */
void addHashSn(Read* r, HashAln** table,
    uint32_t hashSize, bool dupsVerb) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t pos = (a->strand ? a->pos[0] : a->pos[1]);
    uint32_t idx = jenkins_hash_aln(a->chrom, NULL,
      pos, 0, a->strand, 0, SINGLE, hashSize);
    addToHash(a->chrom, NULL, pos, 0, a->strand, 0,
      table, idx, dupsVerb ? r->name : NULL);
  }
}

/* bool checkHashSn()
 * Check a set of singleton alignments for a match in the
 *   hashtable. Return true if *any* match.
 */
bool checkHashSn(Read* r, HashAln** table,
    uint32_t hashSize, File dups, bool dupsVerb,
    bool gzOut) {
  for (int k = 0; k < r->alnLen; k++) {
    Aln* a = r->aln + k;
    uint32_t pos = (a->strand ? a->pos[0] : a->pos[1]);
    uint32_t idx = jenkins_hash_aln(a->chrom, NULL,
      pos, 0, a->strand, 0, SINGLE, hashSize);
    HashAln* h = checkHash(a->chrom, NULL, pos, 0,
      a->strand, 0, SINGLE, table, idx);
    if (h) {
      if (dupsVerb)
        logDup(dups, gzOut, r->name, a->chrom, NULL,
          pos, 0, a->strand, 0, h->name, SINGLE);
      return true;
    }
  }
  return false;
}

/* void findDupsSn()
 * Find PCR duplicates among singleton alignment sets.
 *   Note: the hashtable was already created and
 *   populated by paired and discordant alns.
 */
void findDupsSn(Read** readSn, int readIdxSn,
    int readLenSn, int* countSn, int* dupsSn,
    int* singlePr, bool extendOpt, int extend,
    float asDiff, bool atacOpt, int atacLen5, int atacLen3,
    File bed, bool bedOpt, bool gzOut, bool ctrl,
    int sample, HashAln** table, uint32_t hashSize,
    File dups, bool dupsVerb, uint32_t* order,
    uint32_t* order0, uint32_t* order1, uint32_t* order2,
    uint16_t* qual, uint16_t* qual0, uint16_t* qual1,
    uint16_t* qual2, bool verbose) {

  // sort reads by qual score sum
  uint32_t count = readIdxSn * MAX_SIZE + readLenSn;
  sortReads(readSn, count, order, qual, order0, order1,
    order2, qual0, qual1, qual2);

  // loop through singleton reads
  for (uint32_t i = 0; i < count; i++) {
    Read* r = readSn[order[i] / MAX_SIZE]
      + order[i] % MAX_SIZE;

    // check hashtable for matches
    if (checkHashSn(r, table, hashSize, dups,
        dupsVerb, gzOut))
      (*dupsSn)++;
    else {
      // add alignments to hashtable
      addHashSn(r, table, hashSize, dupsVerb);
      // process alignments too
      *singlePr += processSingle(r->name, r->aln,
        r->alnLen, extendOpt, extend, false,
        NULL, NULL, NULL, NULL,
        r->score, asDiff, r->first,
        atacOpt, atacLen5, atacLen3, bed, bedOpt,
        gzOut, ctrl, sample, verbose);
    }

    // free Read
    free(r->name);
    free(r->aln);

    (*countSn)++;
  }

  // free nodes of hashtable
  for (uint32_t i = 0; i < hashSize; i++) {
    HashAln* tmp;
    HashAln* h = table[i];
    while (h != NULL) {
      free(h->name);
      tmp = h->next;
      free(h);
      h = tmp;
    }
  }
}

/* void findDups()
 * Control elucidation of PCR duplicates. Process reads
 *   that are determined not to be duplicates.
 */
void findDups(Read** readPr, int readIdxPr, int readLenPr,
    Read** readDc, int readIdxDc, int readLenDc,
    Read** readSn, int readIdxSn, int readLenSn,
    HashAln*** table, uint32_t* tableMem,
    HashAln*** tableSn, uint32_t* tableSnMem,
    uint32_t** order, uint32_t** order0, uint32_t** order1,
    uint32_t** order2, uint16_t** qual, uint16_t** qual0,
    uint16_t** qual1, uint16_t** qual2, uint32_t* arrMem,
    int* countPr, int* dupsPr, int* countDc, int* dupsDc,
    int* countSn, int* dupsSn, bool singleOpt,
    int* pairedPr, int* singlePr, double* totalLen,
    bool extendOpt, int extend, bool avgExtOpt,
    float asDiff, bool atacOpt, int atacLen5, int atacLen3,
    File bed, bool bedOpt, bool gzOut, File dups,
    bool dupsVerb, bool ctrl, int sample, bool verbose) {

  // initialize hash table for singletons
  uint32_t hashSizeSn = 0;
  if (singleOpt && (readIdxSn || readLenSn)) {
    // calculate hashtable size
    uint32_t sum = MIN(UINT32_MAX,
      2 * (readIdxPr * MAX_SIZE + readLenPr)
      + 2 * (readIdxDc * MAX_SIZE + readLenDc)
      + readIdxSn * MAX_SIZE + readLenSn);
    hashSizeSn = calcHashSize(sum);

    // create/expand hashtable
    if (hashSizeSn > *tableSnMem) {
      *tableSn = (HashAln**) memrealloc(*tableSn,
        hashSizeSn * sizeof(HashAln*));
      *tableSnMem = hashSizeSn;
    } else
      hashSizeSn = *tableSnMem;
    for (uint32_t i = 0; i < hashSizeSn; i++)
      (*tableSn)[i] = NULL;
  }

  // initialize arrays for efficient sorting
  uint32_t count = MIN(UINT32_MAX,
    MAX(readIdxPr * MAX_SIZE + readLenPr,
    MAX(readIdxDc * MAX_SIZE + readLenDc,
    readIdxSn * MAX_SIZE + readLenSn)));
  if (count > *arrMem) {
    *order = (uint32_t*) memrealloc(*order, count * sizeof(uint32_t));
    *order0 = (uint32_t*) memrealloc(*order0, count * sizeof(uint32_t));
    *order1 = (uint32_t*) memrealloc(*order1, count * sizeof(uint32_t));
    *order2 = (uint32_t*) memrealloc(*order2, count * sizeof(uint32_t));
    *qual = (uint16_t*) memrealloc(*qual, count * sizeof(uint16_t));
    *qual0 = (uint16_t*) memrealloc(*qual0, count * sizeof(uint16_t));
    *qual1 = (uint16_t*) memrealloc(*qual1, count * sizeof(uint16_t));
    *qual2 = (uint16_t*) memrealloc(*qual2, count * sizeof(uint16_t));
    *arrMem = count;
  }

  // evaluate and process paired alignments
  if (readIdxPr || readLenPr)
    findDupsPr(readPr, readIdxPr, readLenPr, countPr,
      dupsPr, pairedPr, totalLen, asDiff, atacOpt,
      atacLen5, atacLen3, bed, bedOpt, gzOut, ctrl,
      sample, table, tableMem, *tableSn, hashSizeSn,
      dups, dupsVerb, *order, *order0, *order1, *order2,
      *qual, *qual0, *qual1, *qual2, verbose);

  if (singleOpt) {
    // with avgExtOpt, calculate average fragment length
    //   and save it as 'extend' with extendOpt=true
    if (avgExtOpt) {
      extend = calcAvgLen(*totalLen, *pairedPr, verbose);
      if (extend)
        extendOpt = true;
    }

    // evaluate and process discordant alignments
    if (readIdxDc || readLenDc)
      findDupsDc(readDc, readIdxDc, readLenDc, countDc,
        dupsDc, singlePr, extendOpt, extend, asDiff,
        atacOpt, atacLen5, atacLen3, bed, bedOpt, gzOut,
        ctrl, sample, table, tableMem, *tableSn,
        hashSizeSn, dups, dupsVerb, *order, *order0,
        *order1, *order2, *qual, *qual0, *qual1, *qual2,
        verbose);

    // evaluate and process singleton alignments
    if (readIdxSn || readLenSn)
      findDupsSn(readSn, readIdxSn, readLenSn, countSn,
        dupsSn, singlePr, extendOpt, extend, asDiff,
        atacOpt, atacLen5, atacLen3, bed, bedOpt, gzOut,
        ctrl, sample, *tableSn, hashSizeSn, dups, dupsVerb,
        *order, *order0, *order1, *order2, *qual, *qual0,
        *qual1, *qual2, verbose);
  }

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
    a->first = true;
  } else {
    a->pos[0] = pnext;
    a->pos[1] = flag & 0x10 ? pos + length : pos;
    a->first = false;
  }

  (*alnLen)++;
  return true;
}

/* bool saveSingleAln()
 * Save the information for an unpaired alignment.
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

/* uint16_t sumQual()
 * Sum an array/string of quality scores.
 */
uint16_t sumQual(char* qual, int len, int offset) {
  if (qual[0] == 0xFF)  // BAM 'null' value
    return 0;
  int sum = 0;
  for (int i = 0; i < len; i++)
    sum += qual[i] - offset;
  return sum > UINT16_MAX ? UINT16_MAX : (uint16_t) sum;
}

/* bool parseAlign()
 * Parse a SAM/BAM alignment record. Save alignment
 *   info to Aln* array. Return true unless the max.
 *   number of alignments has been reached.
 */
bool parseAlign(Aln** aln, int* alnLen, uint16_t flag,
    Chrom* chrom, uint32_t pos, int length, uint32_t pnext,
    int* paired, int* single, int* secPair, int* secSingle,
    int* skipped, bool singleOpt, float score,
    bool dupsOpt, char* qual, int qualLen, int offset,
    uint16_t* qualR1, uint16_t* qualR2) {

  // check for linear template or missing index
  if (flag & 0x1) {
    if ((flag & 0xC0) == 0xC0)
      exit(error("", ERRLINEAR));
    if (!(flag & 0xC0))
      exit(error("", ERRINDEX));
  }

  // save sum of quality scores (only if removing dups)
  if (dupsOpt) {
    if (flag & 0x40) {
      if (! *qualR1 && strcmp(qual, "*"))
        *qualR1 = sumQual(qual, qualLen, offset);
    } else {
      if (! *qualR2 && strcmp(qual, "*"))
        *qualR2 = sumQual(qual, qualLen, offset);
    }
  }

  // paired alignment: save alignment information
  if ((flag & 0x3) == 0x3) {

    // update counts
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
          && (flag & 0x40 ? (! a->first && a->pos[0] == pos)
            : (a->first && a->pos[1] == pos) )
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

  // save to list
  *chrom = (Chrom*) memrealloc(*chrom,
    (*chromLen + 1) * sizeof(Chrom));
  Chrom* c = *chrom + *chromLen;
  c->name = (char*) memalloc(1 + strlen(name));
  strcpy(c->name, name);
  c->len = len;
  c->skip = checkChrom(c->name, xcount, xchrList);
  c->save = ! ctrl; // do not save if ref in ctrl sample only
  c->diff = NULL;
  c->expt = (Pileup*) memalloc(sizeof(Pileup));
  c->expt->end = NULL;
  c->expt->cov = NULL;
  c->exptLen = 0;
  c->exptMem = 0;
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
  c->bed = NULL;
  c->bedLen = 0;
  if (! c->skip)
    saveXBed(c->name, c->len, &c->bedLen, &c->bed,
      xBedLen, xBed, verbose);

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

    // sort order must be queryname
    if (order == NULL || strcmp(order, "queryname"))
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
    bool dupsOpt, File dups, bool dupsVerb, Read*** readPr,
    int* readMemPr, Read*** readDc, int* readMemDc,
    Read*** readSn, int* readMemSn, HashAln*** table,
    uint32_t* tableMem, HashAln*** tableSn,
    uint32_t* tableSnMem, uint32_t** order,
    uint32_t** order0, uint32_t** order1,
    uint32_t** order2, uint16_t** qualA, uint16_t** qual0,
    uint16_t** qual1, uint16_t** qual2, uint32_t* arrMem,
    int* countPr, int* dupsPr, int* countDc, int* dupsDc,
    int* countSn, int* dupsSn, bool verbose) {

  // SAM fields to save
  char* qname, *rname, *cigar, *rnext, *seq, *qual, *extra;
  uint16_t flag;
  uint32_t pos, pnext;
  int32_t tlen;
  uint8_t mapq;

  int alnLen = 0;     // number of alignments for this read
  int unpairIdx = 0;  // \ indexes into unpaired array(s)
  int unpairLen = 0;  // /   (with avgExtOpt)
  int readIdxPr = 0;  // \ indexes into read array readPr
  int readLenPr = 0;  // /   (with dupsOpt)
  int readIdxDc = 0;  // \ indexes into read array readDc
  int readLenDc = 0;  // /   (with dupsOpt)
  int readIdxSn = 0;  // \ indexes into read array readSn
  int readLenSn = 0;  // /   (with dupsOpt)
  uint16_t qualR1 = 0, qualR2 = 0; // sums of quality scores
  bool pastHeader = false;    // to check for misplaced header lines
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
          gzOut, ctrl, sample, dupsOpt,
          readPr, &readIdxPr, &readLenPr, readMemPr,
          readDc, &readIdxDc, &readLenDc, readMemDc,
          readSn, &readIdxSn, &readLenSn, readMemSn,
          qualR1, qualR2, verbose);
      alnLen = 0;
      qualR1 = qualR2 = 0;
      strncpy(readName, qname, MAX_ALNS);
    }

    // save alignment information
    int length = calcDist(qname, seq, cigar); // distance to 3' end
    float score = getScore(extra);
    if (! parseAlign(aln, &alnLen, flag, ref, pos, length,
        pnext, paired, single, secPair, secSingle, skipped,
        singleOpt, score, dupsOpt, qual, strlen(qual),
        SAMQUAL, &qualR1, &qualR2) && verbose)
      fprintf(stderr, "Warning! Read %s has more than %d alignments\n",
        qname, MAX_ALNS);
    // NOTE: the following SAM fields are ignored:
    //   rnext, tlen
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt, unpair,
      &unpairIdx, &unpairLen, unpairMem, asDiff,
      atacOpt, atacLen5, atacLen3, bed, bedOpt,
      gzOut, ctrl, sample, dupsOpt,
      readPr, &readIdxPr, &readLenPr, readMemPr,
      readDc, &readIdxDc, &readLenDc, readMemDc,
      readSn, &readIdxSn, &readLenSn, readMemSn,
      qualR1, qualR2, verbose);

  if (dupsOpt)
    // remove duplicates and process all alignments
    findDups(*readPr, readIdxPr, readLenPr, *readDc,
      readIdxDc, readLenDc, *readSn, readIdxSn, readLenSn,
      table, tableMem, tableSn, tableSnMem, order, order0,
      order1, order2, qualA, qual0, qual1, qual2, arrMem,
      countPr, dupsPr, countDc, dupsDc, countSn, dupsSn,
      singleOpt, pairedPr, singlePr, totalLen, extendOpt,
      extend, avgExtOpt, asDiff, atacOpt, atacLen5,
      atacLen3, bed, bedOpt, gzOut, dups, dupsVerb, ctrl,
      sample, verbose);

  else if (avgExtOpt)
    // process single alignments w/ avgExtOpt
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
    bool ctrl, int sample, bool dupsOpt, File dups,
    bool dupsVerb, Read*** readPr, int* readMemPr,
    Read*** readDc, int* readMemDc, Read*** readSn,
    int* readMemSn, HashAln*** table, uint32_t* tableMem,
    HashAln*** tableSn, uint32_t* tableSnMem,
    uint32_t** order, uint32_t** order0, uint32_t** order1,
    uint32_t** order2, uint16_t** qualA, uint16_t** qual0,
    uint16_t** qual1, uint16_t** qual2, uint32_t* arrMem,
    int* countPr, int* dupsPr, int* countDc, int* dupsDc,
    int* countSn, int* dupsSn, bool verbose) {

  // BAM fields to save
  int32_t refID, pos, l_seq, next_refID, next_pos, tlen;
  uint16_t n_cigar_op, flag;
  uint8_t mapq;
  uint32_t* cigar;
  uint8_t* seq;
  char* read_name, *qual, *extra;

  int alnLen = 0;     // number of alignments for this read
  int unpairIdx = 0;  // \ indexes into unpaired array(s)
  int unpairLen = 0;  // /   (with avgExtOpt)
  int readIdxPr = 0;  // \ indexes into read array readPr
  int readLenPr = 0;  // /   (with dupsOpt)
  int readIdxDc = 0;  // \ indexes into read array readDc
  int readLenDc = 0;  // /   (with dupsOpt)
  int readIdxSn = 0;  // \ indexes into read array readSn
  int readLenSn = 0;  // /   (with dupsOpt)
  uint16_t qualR1 = 0, qualR2 = 0; // sums of quality scores
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
          gzOut, ctrl, sample, dupsOpt,
          readPr, &readIdxPr, &readLenPr, readMemPr,
          readDc, &readIdxDc, &readLenDc, readMemDc,
          readSn, &readIdxSn, &readLenSn, readMemSn,
          qualR1, qualR2, verbose);
      alnLen = 0;
      qualR1 = qualR2 = 0;
      strncpy(readName, read_name, MAX_ALNS);
    }

    // save alignment information
    int length = calcDistBAM(l_seq, n_cigar_op, cigar); // distance to 3' end
    float score = getBAMscore(extra, block_size
      - (int) (extra - line));
    if (! parseAlign(aln, &alnLen, flag, chrom + idx[refID],
        pos, length, next_pos, paired, single, secPair,
        secSingle, skipped, singleOpt, score, dupsOpt,
        qual, l_seq, 0, &qualR1, &qualR2) && verbose)
      fprintf(stderr, "Warning! Read %s has more than %d alignments\n",
        read_name, MAX_ALNS);
    // NOTE: the following BAM fields are ignored:
    //   next_refID, tlen, seq
  }

  // process last set of alns
  if (readName[0] != '\0')
    processAlns(readName, *aln, alnLen, totalLen,
      pairedPr, singlePr, orphan, singleOpt,
      extendOpt, extend, avgExtOpt, unpair,
      &unpairIdx, &unpairLen, unpairMem, asDiff,
      atacOpt, atacLen5, atacLen3, bed, bedOpt,
      gzOut, ctrl, sample, dupsOpt,
      readPr, &readIdxPr, &readLenPr, readMemPr,
      readDc, &readIdxDc, &readLenDc, readMemDc,
      readSn, &readIdxSn, &readLenSn, readMemSn,
      qualR1, qualR2, verbose);

  if (dupsOpt)
    // remove duplicates and process all alignments
    findDups(*readPr, readIdxPr, readLenPr, *readDc,
      readIdxDc, readLenDc, *readSn, readIdxSn, readLenSn,
      table, tableMem, tableSn, tableSnMem, order, order0,
      order1, order2, qualA, qual0, qual1, qual2, arrMem,
      countPr, dupsPr, countDc, dupsDc, countSn, dupsSn,
      singleOpt, pairedPr, singlePr, totalLen, extendOpt,
      extend, avgExtOpt, asDiff, atacOpt, atacLen5,
      atacLen3, bed, bedOpt, gzOut, dups, dupsVerb, ctrl,
      sample, verbose);

  else if (avgExtOpt)
    // process single alignments w/ avgExtOpt
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
    bool dupsOpt, File dups, bool dupsVerb, Read*** readPr,
    int* readMemPr, Read*** readDc, int* readMemDc,
    Read*** readSn, int* readMemSn, HashAln*** table,
    uint32_t* tableMem, HashAln*** tableSn,
    uint32_t* tableSnMem, uint32_t** order,
    uint32_t** order0, uint32_t** order1,
    uint32_t** order2, uint16_t** qual, uint16_t** qual0,
    uint16_t** qual1, uint16_t** qual2, uint32_t* arrMem,
    int* countPr, int* dupsPr, int* countDc, int* dupsDc,
    int* countSn, int* dupsSn, bool verbose) {

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
  char* sortOrder = NULL;
  char* field = strtok(NULL, TAB);
  while (field != NULL) {
    if (!strncmp(field, "SO:", 3))
      sortOrder = field + 3;
    field = strtok(NULL, TAB);
  }
  if (sortOrder == NULL || strcmp(sortOrder, "queryname"))
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

  return parseBAM(in, line, aln, readName, *chromLen,
    *chrom, n_ref, idx, totalLen, unmapped, paired, single,
    pairedPr, singlePr, supp, skipped, lowMapQ, minMapQ,
    secPair, secSingle, orphan, singleOpt, extendOpt,
    extend, avgExtOpt, unpair, unpairMem, asDiff, atacOpt,
    atacLen5, atacLen3, bed, bedOpt, gzOut, ctrl, sample,
    dupsOpt, dups, dupsVerb, readPr, readMemPr, readDc,
    readMemDc, readSn, readMemSn, table, tableMem, tableSn,
    tableSnMem, order, order0, order1, order2, qual, qual0,
    qual1, qual2, arrMem, countPr, dupsPr, countDc, dupsDc,
    countSn, dupsSn, verbose);
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
 * Load genomic regions to exclude from BED file(s).
 *   Return number saved.
 */
int loadBED(char* xFile, char* line, Bed** xBed) {

  // loop through BED files
  int count = 0;  // count of intervals saved
  char* end;
  char* filename = strtok_r(xFile, COM, &end);
  while (filename) {

    // open BED file
    File in;
    bool gz = openRead(filename, &in);

    // load BED records
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
        char msg[MAX_ALNS];
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
      exit(error(filename, ERRCLOSE));

    filename = strtok_r(NULL, COM, &end);
  }

  return count;
}

/* void findPeaksOnly()
 * Control peak-calling directly from logfile (-f).
 */
void findPeaksOnly(char* logFile, char* outFile,
    bool gzOut, int xcount, char** xchrList, char* xFile,
    float pqvalue, bool qvalOpt, int minLen, int maxGap,
    float minAUC, bool verbose) {

  // save genomic regions to exclude
  char* line = (char*) memalloc(MAX_SIZE);
  Bed* xBed = NULL;
  int xBedLen = 0;
  if (xFile != NULL)
    xBedLen = loadBED(xFile, line, &xBed);

  // open files
  File in, out;
  bool gz = openRead(logFile, &in);
  openWrite(outFile, &out, gzOut);
  if (verbose)
    fprintf(stderr, "Peak-calling from log file: %s\n",
      logFile);

  // call peaks
  callPeaksLog(in, gz, out, gzOut, line, xcount, xchrList,
    xBedLen, xBed, pqvalue, qvalOpt, minLen, maxGap,
    minAUC, verbose);

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
  free(line);

  // close files
  if ( ( gz && gzclose(in.gzf) != Z_OK )
      || ( ! gz && fclose(in.f) ) )
    exit(error(logFile, ERRCLOSE));
  if ( ( gzOut && gzclose(out.gzf) != Z_OK )
      || ( ! gzOut && fclose(out.f) ) )
    exit(error(outFile, ERRCLOSE));
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
    int atacLen, bool dupsOpt, int countPr, int dupsPr,
    int countDc, int dupsDc, int countSn, int dupsSn) {
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
  if (dupsOpt) {
    fprintf(stderr, "  PCR duplicates --\n");
    fprintf(stderr, "    Paired aln sets:    %10d\n", countPr);
    fprintf(stderr, "      duplicates:       %10d (%.1f%%)\n",
      dupsPr, countPr ? 100.0f * dupsPr / countPr : 0.0f);
    if (singleOpt) {
      fprintf(stderr, "    Discordant aln sets:%10d\n", countDc);
      fprintf(stderr, "      duplicates:       %10d (%.1f%%)\n",
        dupsDc, countDc ? 100.0f * dupsDc / countDc : 0.0f);
      fprintf(stderr, "    Singleton aln sets: %10d\n", countSn);
      fprintf(stderr, "      duplicates:       %10d (%.1f%%)\n",
        dupsSn, countSn ? 100.0f * dupsSn / countSn : 0.0f);
    }
  }
  fprintf(stderr, "  Fragments analyzed:   %10d\n", singlePr + pairedPr);
  fprintf(stderr, "    Full fragments:     %10d\n", pairedPr);
  if (pairedPr && ! atacOpt)
    fprintf(stderr, "      (avg. length: %.1fbp)\n", avgLen);
  if (singleOpt) {
    fprintf(stderr, "    Half fragments:     %10d\n", singlePr);
    if (singlePr) {
      fprintf(stderr, "      (from unpaired alns");
      if (extendOpt)
        fprintf(stderr, ", extended to %dbp", extend);
      else if (avgExtOpt && pairedPr)
        fprintf(stderr, ", extended to %dbp", (int) (avgLen + 0.5));
      fprintf(stderr, ")\n");
    }
  }
  if (atacOpt) {
    fprintf(stderr, "    ATAC-seq cut sites: %10d\n",
      2 * pairedPr + singlePr);
    fprintf(stderr, "      (expanded to length %dbp)\n", atacLen);
  }
}

/* void runProgram()
 * Controls the opening/closing of files, and parsing
 *   of input files by readSAM() or readBAM().
 *   Pileup values are computed by savePileupExpt() or
 *   savePileupCtrl(), and p-values for each experimental/
 *   control pair are calculated by savePval().
 *   Results for all replicates are passed to findPeaks().
 * If calling peaks only (from logfile), pass control
 *   directly to findPeaksOnly().
 */
void runProgram(char* inFile, char* ctrlFile, char* outFile,
    char* logFile, char* pileFile, char* bedFile,
    bool gzOut, bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt, int minMapQ, int xcount,
    char** xchrList, char* xFile, float pqvalue,
    bool qvalOpt, int minLen, int maxGap, float minAUC,
    float asDiff, bool atacOpt, int atacLen5, int atacLen3,
    bool dupsOpt, char* dupsFile, bool peaksOpt,
    bool peaksOnly, bool verbose) {

  // option to call peaks only, from already produced log file
  if (peaksOnly) {
    findPeaksOnly(logFile, outFile, gzOut, xcount,
      xchrList, xFile, pqvalue, qvalOpt, minLen, maxGap,
      minAUC, verbose);
    return;
  }

  // open optional output files
  File bed, pile, dups;
  if (bedFile != NULL)
    openWrite(bedFile, &bed, gzOut);
  if (pileFile != NULL)
    openWrite(pileFile, &pile, gzOut);
  bool dupsVerb = false;
  if (dupsOpt && dupsFile != NULL) {
    openWrite(dupsFile, &dups, gzOut);
    dupsVerb = true;
  }

  // initialize variables
  char* line = (char*) memalloc(MAX_SIZE);
  int chromLen = 0;     // number of reference sequences
  Chrom* chrom = NULL;  // array of reference sequences
  Aln* aln = (Aln*) memalloc(MAX_ALNS * sizeof(Aln)); // array of saved alns
  int unpairMem = 0;    // number of unpaired alns (for avg-ext option)
  Aln** unpair = NULL;  // array of unpaired alns (for avg-ext option)
  char* readName = memalloc(MAX_ALNS + 1);  // name of read being analyzed
  readName[0] = readName[MAX_ALNS] = '\0';
  int sample = 0;       // number of sample pairs analyzed

  // variables for duplicate removal option
  Read** readPr = NULL;     // array of reads with properly paired aln(s)
  int readMemPr = 0;        // index into readPr array
  Read** readDc = NULL;     // array of reads with discordant aln(s)
  int readMemDc = 0;        // index into readDc array
  Read** readSn = NULL;     // array of reads with singleton aln(s)
  int readMemSn = 0;        // index into readSn array
  HashAln** table = NULL;   // hashtable for paired/discordant alns (recycled)
  uint32_t tableMem = 0;    // size of above hashtable
  HashAln** tableSn = NULL; // hashtable for singletons alns
  uint32_t tableSnMem = 0;  // size of above hashtable
  uint32_t* order = NULL;   // array for order of reads to be processed
  uint32_t* order0 = NULL;  // |
  uint32_t* order1 = NULL;  // | temp arrays for order
  uint32_t* order2 = NULL;  // |
  uint16_t* qual = NULL;    // array of quality score sums (sorting key)
  uint16_t* qual0 = NULL;   // |
  uint16_t* qual1 = NULL;   // | temp arrays for qual
  uint16_t* qual2 = NULL;   // |
  uint32_t arrMem = 0;      // length of above order/qual arrays

  // save genomic regions to exclude
  Bed* xBed = NULL;
  int xBedLen = 0;
  if (xFile != NULL)
    xBedLen = loadBED(xFile, line, &xBed);

  // loop through input files (experimental and control)
  char* end1, *end2;
  char* exptName = strtok_r(inFile, COM, &end1);
  char* ctrlName = ctrlFile == NULL ? NULL
    : strtok_r(ctrlFile, COM, &end2);
  while (exptName) {

    // reset 'save' bools of each Chrom
    for (int j = 0; j < chromLen; j++)
      (chrom + j)->save = false;

    // process matching experimental/control files
    double fragLen = 0.0; // total weighted length of all experimental fragments
    for (int i = 0; i < 2; i++) {

      // get expt/ctrl filename
      char* filename = exptName;
      if (i) {
        filename = ctrlName;
        if (ctrlName != NULL && !strcmp(ctrlName, "null"))
          filename = NULL;
        if (filename == NULL) {
          if (verbose)
            fprintf(stderr, "- control file #%d not provided -\n",
              sample);
          savePileupNoCtrl(chrom, chromLen, fragLen,
            verbose);
          break;
        }
      }

      // open input file
      File in;
      bool gz = openRead(filename, &in);
      bool bam = checkBAM(in, gz);
      if (verbose)
        fprintf(stderr, "Processing %s file #%d: %s\n",
          i ? "control" : "experimental", sample, filename);
      if (dupsVerb) {
        if (gzOut)
          gzprintf(dups.gzf, "# %s file #%d: %s\n",
            i ? "control" : "experimental", sample, filename);
        else
          fprintf(dups.f, "# %s file #%d: %s\n",
            i ? "control" : "experimental", sample, filename);
      }

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
        lowMapQ = 0, secPair = 0, secSingle = 0,
        countPr = 0, dupsPr = 0, countDc = 0, dupsDc = 0,
        countSn = 0, dupsSn = 0;  // counting variables
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
          bed, bedFile != NULL, gzOut, i, sample, dupsOpt,
          dups, dupsVerb, &readPr, &readMemPr, &readDc,
          &readMemDc, &readSn, &readMemSn, &table,
          &tableMem, &tableSn, &tableSnMem, &order,
          &order0, &order1, &order2, &qual, &qual0, &qual1,
          &qual2, &arrMem, &countPr, &dupsPr, &countDc,
          &dupsDc, &countSn, &dupsSn, verbose);
      else
        count = readSAM(in, gz, line, &aln, readName,
          &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ,
          minMapQ, xcount, xchrList, xBedLen, xBed,
          &secPair, &secSingle, &orphan, &chromLen, &chrom,
          singleOpt, extendOpt, extend, avgExtOpt, &unpair,
          &unpairMem, asDiff, atacOpt, atacLen5, atacLen3,
          bed, bedFile != NULL, gzOut, i, sample, dupsOpt,
          dups, dupsVerb, &readPr, &readMemPr, &readDc,
          &readMemDc, &readSn, &readMemSn, &table,
          &tableMem, &tableSn, &tableSnMem, &order,
          &order0, &order1, &order2, &qual, &qual0, &qual1,
          &qual2, &arrMem, &countPr, &dupsPr, &countDc,
          &dupsDc, &countSn, &dupsSn, verbose);

      // log counts
      if (verbose)
        logCounts(count, unmapped, supp, skipped, chrom,
          chromLen, minMapQ, lowMapQ, paired, secPair,
          orphan, single, secSingle, singlePr, pairedPr,
          totalLen, singleOpt, extendOpt, extend,
          avgExtOpt, bam, atacOpt, atacLen5 + atacLen3,
          dupsOpt, countPr, dupsPr, countDc, dupsDc,
          countSn, dupsSn);

      // save pileup values
      if (i)
        savePileupCtrl(chrom, chromLen, fragLen, verbose);
      else
        fragLen = savePileupExpt(chrom, chromLen);

      // close input file
      if ( (gz && gzclose(in.gzf) != Z_OK)
          || (! gz && fclose(in.f)) )
        exit(error(filename, ERRCLOSE));
    }

    // calculate p-values
    savePval(chrom, chromLen, sample, pile,
      pileFile != NULL, exptName, ctrlName, gzOut);

    exptName = strtok_r(NULL, COM, &end1);
    ctrlName = ctrlFile == NULL ? NULL
      : strtok_r(NULL, COM, &end2);
    sample++;
  }

  // free 'diff' arrays
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (chr->diff) {
      free(chr->diff->frac);
      free(chr->diff->cov);
      free(chr->diff);
    }
  }

  // open output files
  File out, log;
  if (peaksOpt)
    openWrite(outFile, &out, gzOut);
  if (logFile != NULL)
    openWrite(logFile, &log, gzOut);

  // find peaks
  findPeaks(out, log, logFile != NULL, gzOut, chrom,
    chromLen, &sample, pqvalue, qvalOpt, minLen,
    maxGap, minAUC, peaksOpt, verbose);

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
  if (dupsOpt) {
    for (int i = 0; i < readMemPr; i++)
      free(readPr[i]);
    free(readPr);
    for (int i = 0; i < readMemDc; i++)
      free(readDc[i]);
    free(readDc);
    for (int i = 0; i < readMemSn; i++)
      free(readSn[i]);
    free(readSn);
    free(table);
    free(tableSn);
    free(order);
    free(order0);
    free(order1);
    free(order2);
    free(qual);
    free(qual0);
    free(qual1);
    free(qual2);
  }
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom + i;
    if (! chr->skip) {
      if (chr->bedLen)
        free(chr->bed);
      if (qvalOpt && chr->qval) {
        free(chr->qval->end);
        free(chr->qval->cov);
        free(chr->qval);
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
      free(chr->expt->end);
      free(chr->expt->cov);
    }
    free(chr->ctrl);
    free(chr->expt);
    free(chr->name);
  }
  free(chrom);
  free(aln);
  free(readName);
  free(line);

  // close files
  if (peaksOpt && ( ( gzOut && gzclose(out.gzf) != Z_OK )
      || ( ! gzOut && fclose(out.f) ) ) )
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
  if (dupsVerb && ( ( ! gzOut && fclose(dups.f) )
      || ( gzOut && gzclose(dups.gzf) != Z_OK ) ) )
    exit(error(dupsFile, ERRCLOSE));
}

/* int saveXChrom()
 * Save list of chromosomes (ref names) to exclude.
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
    *xFile = NULL, *dupsFile = NULL;
  char* xchrom = NULL;
  int extend = 0, minMapQ = 0, minLen = DEFMINLEN,
    maxGap = DEFMAXGAP, atacLen5 = DEFATAC, atacLen3 = 0;
  float asDiff = 0.0f, pqvalue = DEFQVAL, minAUC = DEFAUC;
  bool singleOpt = false, extendOpt = false,
    avgExtOpt = false, atacOpt = false, gzOut = false,
    qvalOpt = true, dupsOpt = false,
    peaksOpt = true, peaksOnly = false;
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
      case UNPAIROPT: singleOpt = true; break;
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
      case MINAUC: minAUC = getFloat(optarg); break;
      case MINLEN: minLen = getInt(optarg); break;
      case MAXGAP: maxGap = getInt(optarg); break;
      case DUPSOPT: dupsOpt = true; break;
      case DUPSFILE: dupsFile = optarg; break;
      case NOPEAKS: peaksOpt = false; break;
      case PEAKSONLY: peaksOnly = true; break;
      case VERBOSE: verbose = true; break;
      case VERSOPT: printVersion(); break;
      case HELP: usage(); break;
      default: exit(EXIT_FAILURE);
    }
  if (optind < argc)
    exit(error(argv[optind], ERRPARAM));

  // check for argument errors
  if ((peaksOpt && outFile == NULL)
      || (peaksOnly && logFile == NULL)
      || (!peaksOnly && inFile == NULL)) {
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
    avgExtOpt = extendOpt = false;  // no unpaired extensions in ATAC-seq mode
    if (atacLen5 <= 0)
      exit(error("", ERRATAC));
    // split atacLen into atacLen5 and atacLen3
    atacLen3 = (int) (atacLen5 / 2.0f + 0.5f);  // round up for 3' end
    atacLen5 /= 2;
  }
  if (minLen < 0)
    exit(error("", ERRMINLEN));
  if (minAUC < 0.0f)
    exit(error("", ERRMINAUC));
  if (asDiff < 0.0f)
    exit(error("", ERRASDIFF));

  // save list of chromosomes to exclude
  int xcount = 0;
  char** xchrList = NULL;
  if (xchrom != NULL)
    xcount = saveXChrom(xchrom, &xchrList);

  // adjust significance level to -log scale
  if (pqvalue <= 0.0f || pqvalue > 1.0f)
    exit(error("", ERRPQVAL));
  pqvalue = -log10f(pqvalue);

  // send arguments to runProgram()
  runProgram(inFile, ctrlFile, outFile, logFile, pileFile,
    bedFile, gzOut, singleOpt, extendOpt, extend,
    avgExtOpt, minMapQ, xcount, xchrList, xFile, pqvalue,
    qvalOpt, minLen, maxGap, minAUC, asDiff,
    atacOpt, atacLen5, atacLen3, dupsOpt, dupsFile,
    peaksOpt, peaksOnly, verbose);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getArgs(argc, argv);
  return EXIT_SUCCESS;
}
