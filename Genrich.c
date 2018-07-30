/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.0
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
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
  fprintf(stderr, "  -%c  <file>       Output BED file\n", OUTFILE);
  fprintf(stderr, "Optional arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input SAM/BAM file(s) for control sample(s)\n", CTRLFILE);
  fprintf(stderr, "Filtering options:\n");
  fprintf(stderr, "  -%c  <arg>        Comma-separated list of chromosomes to ignore\n", XCHROM);
  fprintf(stderr, "  -%c  <int>        Minimum MAPQ to keep an alignment (def. 0)\n", MINMAPQ);
  fprintf(stderr, "Options for unpaired alignments:\n");
  fprintf(stderr, "  -%c               Print unpaired alignments (def. false)\n", SINGLEOPT);
  fprintf(stderr, "  -%c  <int>        Print unpaired alignments, with fragment length\n", EXTENDOPT);
  fprintf(stderr, "                     increased to specified value\n");
  fprintf(stderr, "  -%c               Print unpaired alignments, with fragment length\n", AVGEXTOPT);
  fprintf(stderr, "                     increased to average value of paired alignments\n");
  fprintf(stderr, "I/O options:\n");
  fprintf(stderr, "  -%c               Option to gzip output(s)\n", GZOPT);
/*  fprintf(stderr, "  -%c               Option to check for dovetailing (with 3' overhangs)\n", DOVEOPT);
  fprintf(stderr, "  -%c  <int>        Minimum overlap of dovetailed alignments (def. %d)\n", DOVEOVER, DEFDOVE);
  fprintf(stderr, "  -%c  <file>       Log file for stitching results of each read pair\n", LOGFILE);
  fprintf(stderr, "  -%c  <file>       FASTQ files for reads that failed stitching\n", UNFILE);
  fprintf(stderr, "                     (output as <file>%s and <file>%s)\n", ONEEXT, TWOEXT);
  fprintf(stderr, "  -%c  <file>       Log file for dovetailed reads (adapter sequences)\n", DOVEFILE);
  fprintf(stderr, "  -%c  <file>       Log file for formatted alignments of merged reads\n", ALNFILE);
  fprintf(stderr, "  -%c               Option to produce interleaved FASTQ output(s)\n", INTEROPT);
  fprintf(stderr, "  -%c  <file>       Use given error profile for merged qual scores\n", QUALFILE);
  fprintf(stderr, "  -%c               Use 'fastq-join' method for merged qual scores\n", FJOINOPT);
  fprintf(stderr, "  -%c  <int>        FASTQ quality offset (def. %d)\n", QUALITY, OFFSET);
  fprintf(stderr, "  -%c  <int>        Maximum input quality score (0-based; def. %d)\n", SETQUAL, MAXQUAL);
  fprintf(stderr, "  -%c  <int>        Number of threads to use (def. %d)\n", THREADS, DEFTHR);
  fprintf(stderr, "  -%c               Option to print status updates/counts to stderr\n", VERBOSE);
*/
  exit(-1);
}

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

/*** BED printing ***/

/* void printInterval()
 * Print BED interval, with pileup as the value.
 */
void printInterval(File out, bool gzOut, Chrom* chrom,
    int start, int end, float treatVal, float ctrlVal) {
  if (gzOut)
    gzprintf(out.gzf, "%s\t%d\t%d\t%.5f\t%.5f\n",
      chrom->name, start, end, treatVal, ctrlVal);
  else
    fprintf(out.f, "%s\t%d\t%d\t%.5f\t%.5f\n",
      chrom->name, start, end, treatVal, ctrlVal);
}

/* void printPileup()
 * Controls printing of pileup values for each Chrom.
 */
void printPileup(File out, bool gzOut, Chrom** chrom,
    int chromLen) {

  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom[i];
    if (chr->skip)
      continue;

    if (chr->treat == NULL) {
      // no treatment coverage -- shortcut printing
      if (chr->ctrl == NULL)
        printInterval(out, gzOut, chr, 0, chr->len, 0, 0);
      else {
        int start = 0;
        for (int j = 0; j < chr->ctrlLen; j++) {
          printInterval(out, gzOut, chr,
            start, chr->ctrl->end[j],
            0, chr->ctrl->cov[j]);
          start = chr->ctrl->end[j];
        }
      }
      continue;
    }

    // print output from treatment and control
    int start = 0;
    int k = 0;
    for (int j = 0; j < chr->treatLen; j++) {
      while (chr->ctrl->end[k] < chr->treat->end[j]) {
        printInterval(out, gzOut, chr,
          start, chr->ctrl->end[k],
          chr->treat->cov[j], chr->ctrl->cov[k]);
        start = chr->ctrl->end[k];
        k++;
      }
      printInterval(out, gzOut, chr,
        start, chr->treat->end[j],
        chr->treat->cov[j], chr->ctrl->cov[k]);
      if (chr->ctrl->end[k] == chr->treat->end[j])
        k++;
      start = chr->treat->end[j];
    }
  }
}

/* float calcLambda()
 * Calculate a background lambda value: sum of fragment
 *   lengths divided by total genome length.
 */
float calcLambda(Chrom** chrom, int chromLen,
    unsigned long fragLen) {
  unsigned long genomeLen = 0;
  for (int i = 0; i < chromLen; i++)
    if (!chrom[i]->skip)
      genomeLen += chrom[i]->len;
  if (! genomeLen)
    exit(error("", ERRGEN));
//printf("genomeLen is %ld; fragLen is %ld; lambda = %.9f\n",
//  genomeLen, fragLen, fragLen / (double) genomeLen);
  return fragLen / (double) genomeLen;
}

/* void savePileupNoCtrl()
 * When no control is available, save the control
 *   pileup as the background lambda value.
 */
void savePileupNoCtrl(Chrom** chrom, int chromLen, long fragLen) {
  float lambda = calcLambda(chrom, chromLen, fragLen);
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom[i];
    if (chr->skip)
      continue;

    chr->ctrl = (Pileup*) memalloc(sizeof(Pileup));
    chr->ctrl->end = (int*) memalloc(sizeof(int));
    chr->ctrl->end[0] = chr->len;
    chr->ctrl->cov = (float*) memalloc(sizeof(float));
    chr->ctrl->cov[0] = lambda;
    chr->ctrlLen = 1;
  }
}

/* unsigned long savePileup()
 * Save pileup values for each Chrom from diff arrays
 *   and background lambda value.
 *   Return total length of all fragments (weighted).
 */
unsigned long savePileup(Chrom** chrom, int chromLen,
    unsigned long fragLen, bool ctrl) {

  // calculate background lambda value
  float lambda = 0.0;
  if (ctrl)
    lambda = calcLambda(chrom, chromLen, fragLen);

  // create pileup for each chrom
  fragLen = 0;
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom[i];
    if (chr->skip)
      continue;
    if (chr->diff == NULL)
      continue;

    // saving control or pileup values?
    Pileup** p = ctrl ? &chr->ctrl : &chr->treat;
    int* size = ctrl ? &chr->ctrlLen : &chr->treatLen;

    // determine number of pileup intervals
    int num = 1;
    for (int i = 1; i < chr->len; i++)
      if (chr->diff[i] != 0)
        num++;

    // create pileup arrays
    *p = (Pileup*) memalloc(sizeof(Pileup));
    (*p)->end = (int*) memalloc(num * sizeof(int));
    (*p)->cov = (float*) memalloc(num * sizeof(float));
    *size = num;

    // save pileup values
    int pos = 0;              // position in pileup arrays
    float val = chr->diff[0]; // pileup value
    int start = 0;            // beginning coordinate of interval
    int i;
    for (i = 1; i < chr->len; i++)
      if (chr->diff[i] != 0) {
        (*p)->end[pos] = i;
        (*p)->cov[pos] = MAX(val, lambda);

        // keep track of fragment length
        fragLen += (i - start) * val;
        start = i;

        // update pileup value
        val += chr->diff[i];
        pos++;
      }

    // save final interval
    (*p)->end[pos] = i;
    (*p)->cov[pos] = MAX(val, lambda);
    fragLen += (i - start) * val;
  }

  return fragLen;
}

/*** Save SAM/BAM header info ***/

/* int saveChrom()
 * If chromosome (reference sequence) has not been
 *   saved yet, save it to the array. Return the index.
 */
int saveChrom(char* name, int len, int* chromLen,
    Chrom*** chrom, int xcount, char** xchrList) {
  // determine if chrom has been saved already
  for (int i = 0; i < *chromLen; i++)
    if (!strcmp((*chrom)[i]->name, name))
      return i;

  // determine if chrom is on skipped list
  bool skip = false;
  for (int i = 0; i < xcount; i++)
    if (!strcmp(xchrList[i], name)) {
      skip = true;
      break;
    }

  // create Chrom*
  Chrom* c = (Chrom*) memalloc(sizeof(Chrom));
  c->name = (char*) memalloc(1 + strlen(name));
  strcpy(c->name, name);
  c->len = len;
  c->skip = skip;
  c->diff = NULL;
  c->treat = NULL;
  c->ctrl = NULL;

  // save to list
  *chrom = (Chrom**) memrealloc(*chrom,
    (*chromLen + 1) * sizeof(Chrom*));
  (*chrom)[*chromLen] = c;
  (*chromLen)++;
  return *chromLen - 1;
}

/*** Alignment parsing / conversion to intervals ***/

/* int saveInterval()
 * Save BED interval for a read/fragment.
 *   Return fragment length.
 */
int saveInterval(Chrom* chrom, int start, int end,
    char* qname) {
  // check validity of positions
  if (start < 0) {
    fprintf(stderr, "Warning! Read %s prevented from extending below 0 on %s\n",
      qname, chrom->name);
    start = 0;
  }
  if (end > chrom->len) {
    fprintf(stderr, "Warning! Read %s prevented from extending past %d on %s\n",
      qname, chrom->len, chrom->name);
    end = chrom->len;
  }

  // save ends of fragment to 'diff' array
  if (chrom->diff == NULL) {
    chrom->diff = (float*) memalloc(chrom->len * sizeof(float));
    for (int i = 0; i < chrom->len; i++)
      chrom->diff[i] = 0.0;
  }
  chrom->diff[start] += 1;
  chrom->diff[end] -= 1;

  return end - start;
}


/* void saveSingle()
 * Control printing for an unpaired alignment.
 */
void saveSingle(Chrom* chrom, char* qname, uint16_t flag,
    uint32_t pos, int length, bool extendOpt, int extend) {
  if (extendOpt) {
    if (flag & 0x10)
      saveInterval(chrom, pos + length - extend, pos + length, qname);
    else
      saveInterval(chrom, pos, pos + extend, qname);
  } else
    saveInterval(chrom, pos, pos + length, qname);
}

/* int saveAvgExt()
 * Save complete intervals for unpaired alignments
 *   with "extend to average length" option.
 *   Return number printed.
 */
int saveAvgExt(int readLen, Read** unpaired,
    unsigned long totalLen, int pairedPr) {
  // determine average fragment length
  int avgLen = 0;
  if (! pairedPr) {
    fprintf(stderr, "Warning! No paired alignments to calculate avg ");
    fprintf(stderr, "frag length --\n  Printing singletons \"as is\"\n");
  } else
    avgLen = (int) ((double) totalLen / pairedPr + 0.5);

  int printed = 0;  // counting variable
  for (int i = 0; i < readLen; i++) {
    Read* r = unpaired[i];
    if (! avgLen)
      saveInterval(r->chrom, r->pos[0], r->pos[1], r->name);
    else if (r->strand)
      saveInterval(r->chrom, r->pos[0], r->pos[0] + avgLen, r->name);
    else
      saveInterval(r->chrom, r->pos[1] - avgLen, r->pos[1], r->name);
    printed++;

    // free memory
    free(r->name);
    free(r);
  }
  free(unpaired);

  return printed;
}

/* void saveAvgExtList()
 * Save info for an unpaired alignment to list
 *   (for "extend to average length" option).
 */
void saveAvgExtList(int* readLen, int* readMem,
    Read*** unpaired, char* qname, uint16_t flag,
    Chrom* chrom, uint32_t pos, int length) {
  // create new Read
  Read* r = (Read*) memalloc(sizeof(Read));
  r->name = (char*) memalloc(1 + strlen(qname));
  strcpy(r->name, qname);
  r->chrom = chrom;
  r->strand = flag & 0x10 ? false : true;
  r->pos[0] = pos;
  r->pos[1] = pos + length;

  // save to list
  if (*readLen + 1 > *readMem) {
    *unpaired = (Read**) memrealloc(*unpaired,
      (*readMem + 10000) * sizeof(Read*));
    *readMem += 10000;
  }
  (*unpaired)[*readLen] = r;
  (*readLen)++;
}

/* int saveFragment()
 * Save full fragment for a proper pair. Return length.
 */
int saveFragment(Chrom* chrom, Read* r) {
  // ensure start < end
  int start, end;
  if (r->pos[0] > r->pos[1]) {
    start = r->pos[1];
    end = r->pos[0];
  } else {
    start = r->pos[0];
    end = r->pos[1];
  }
  return saveInterval(chrom, start, end, r->name);
}

/* Read* savePaired()
 * Save the position for a properly paired alignment.
 *   If its pair has been analyzed, return the previous
 *   read (for easy removal from the linked list).
 */
Read* savePaired(Read* dummy, char* qname, uint16_t flag,
    Chrom* chrom, uint32_t pos) {
  // determine if read has been analyzed
  Read* r, *prev = dummy;
  for (r = dummy->next; r != NULL; r = r->next) {
    if (! strcmp(qname, r->name))
      break;
    prev = r;
  }

  // if analyzed, save pos and return match
  if (r != NULL) {
    if (flag & 0x40) {
      if (r->pos[0] != -1)
        exit(error(r->name, ERRREP));
      r->pos[0] = pos;
    } else {
      if (r->pos[1] != -1)
        exit(error(r->name, ERRREP));
      r->pos[1] = pos;
    }
    return prev;
  }

  // create new Read
  r = (Read*) memalloc(sizeof(Read));
  r->name = (char*) memalloc(1 + strlen(qname));
  strcpy(r->name, qname);
  r->chrom = chrom;
  r->paired = false;

  // save position
  if (flag & 0x40) {
    r->pos[0] = pos;
    r->pos[1] = -1;
  } else {
    r->pos[0] = -1;
    r->pos[1] = pos;
  }

  // insert r into linked list
  r->next = dummy->next;
  dummy->next = r;
  return NULL;
}

/* void parseAlign()
 * Parse a SAM/BAM alignment record. For properly paired
 *   alignments, save BED interval only if both alignments
 *   have been analyzed. For singleton alignments, save
 *   BED interval only if desired; for avg-ext option, save
 *   alignment until extension length can be calculated from
 *   paired alignments.
 */
void parseAlign(File out, bool gzOut,
    int* readLen, int* readMem, Read*** unpaired, Read* dummy,
    char* qname,
    uint16_t flag, Chrom* chrom, uint32_t pos, int length,
    unsigned long* totalLen, int* paired, int* single, int* pairedPr,
    int* singlePr, bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt) {

  if ((flag & 0x3) == 0x3) {
    // paired alignment: save information
    (*paired)++;
    Read* prev = savePaired(dummy, qname, flag, chrom,
      flag & 0x10 ? pos + length : pos);
    if (prev != NULL) {
      // both alignments analyzed: print BED interval
      Read* r = prev->next;
      (*pairedPr)++;
      *totalLen += saveFragment(chrom, r);
      r->paired = true;

      // remove Read from linked list
      prev->next = r->next;
      free(r->name);
      free(r);
    }
  } else {
    // unpaired alignment
    (*single)++;
    if (singleOpt) {
      if (avgExtOpt)
        // for average-extension option, save alignment
        saveAvgExtList(readLen, readMem, unpaired, qname,
          flag, chrom, pos, length);
      else {
        saveSingle(chrom, qname, flag, pos, length,
          extendOpt, extend);
        (*singlePr)++;
      }
    }
  }
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
      case POS: *pos = getInt(field) - 1; break;
      case MAPQ: *mapq = getInt(field); break;
      case CIGAR: *cigar = field; break;
      case RNEXT: *rnext = field; break;
      case PNEXT: *pnext = getInt(field); break;
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

/* void loadChrom()
 * Check SAM header line for chromosome info.
 */
void loadChrom(char* line, int* chromLen, Chrom*** chrom,
    int xcount, char** xchrList) {
  // determine if SAM header line has chrom info
  char* tag = strtok(line, TAB);
  if (tag == NULL || strcmp(tag, "@SQ"))
    return;
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
    xcount, xchrList);
}

/* int readSAM()
 * Parse the alignments in a SAM file.
 */
int readSAM(File in, bool gz, File out, bool gzOut, char* line,
    File un1, File un2, bool unOpt, File log,
    bool logOpt, int overlap, bool dovetail, int doveOverlap,
    File dove, bool doveOpt, File aln, int alnOpt,
    float mismatch, bool maxLen,
    unsigned long* totalLen, int* unmapped, int* paired, int* single,
    int* pairedPr, int* singlePr,
    int* supp, int* skipped, int* lowMapQ,
    int minMapQ,
    int xcount, char** xchrList, int* chromLen, Chrom*** chrom,
    int* readLen, int* readMem, Read** unpaired, Read* dummy,
    bool singleOpt, bool extendOpt, int extend, bool avgExtOpt,
    int threads) {

  // SAM fields to save
  char* qname, *rname, *cigar, *rnext, *seq, *qual, *extra;
  uint16_t flag;
  uint32_t pos, pnext;
  int32_t tlen;
  uint8_t mapq;

  int count = 0;
  while (getLine(line, MAX_SIZE, in, gz) != NULL) {

    if (line[0] == '@') {
      // load chrom lengths from header
      loadChrom(line, chromLen, chrom, xcount, xchrList);
      continue;
    }

    // parse SAM record
    qname = strtok(line, TAB);
    if (qname == NULL)
      exit(error(line, ERRSAM));
    if (! loadFields(&flag, &rname, &pos, &mapq, &cigar,
        &rnext, &pnext, &tlen, &seq, &qual, &extra))
      exit(error(line, ERRSAM));

    count++;
    if (flag & 0x4) {
      // skip unmapped
      (*unmapped)++;
      continue;
    }
    if (! strcmp(qname, "*") || ! strcmp(rname, "*")
        || pos < 0)
      // insufficient alignment info
      exit(error(line, ERRSAM));
    if (flag & 0xF00) {
      // skip supplementary/secondary alignments
      (*supp)++;
      continue;
    }
    // find matching Chrom (reference sequence)
    Chrom* ref = NULL;
    for (int i = 0; i < *chromLen; i++)
      if (! strcmp((*chrom)[i]->name, rname)) {
        ref = (*chrom)[i];
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

    // save alignment information
    int length = calcDist(qname, seq, cigar); // distance to 3' end
    parseAlign(out, gzOut, readLen, readMem, &unpaired, dummy,
      qname, flag, ref, pos, length, totalLen,
      paired, single, pairedPr, singlePr,
      singleOpt, extendOpt, extend, avgExtOpt);
    // NOTE: the following SAM fields are ignored:
    //   rnext, pnext, tlen, qual, extra (optional fields)
  }

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    *singlePr += saveAvgExt(*readLen, unpaired,
      *totalLen, *pairedPr);


/*******************************************/

/*
  // initialize omp locks -- out, un, log, dove, aln
  omp_lock_t lock[OMP_LOCKS];
  for (int i = 0; i < OMP_LOCKS; i++)
    omp_init_lock(&lock[i]);

  // process files in parallel
  int count = 0, stitchRed = 0;
  #pragma omp parallel num_threads(threads) reduction(+: count, stitchRed)
  {

    // allocate memory for both reads
    char** read1 = (char**) memalloc(FASTQ * sizeof(char*));
    char** read2 = (char**) memalloc((FASTQ + EXTRA) * sizeof(char*));
    for (int i = 0; i < FASTQ + EXTRA; i++) {
      if (i < FASTQ)
        read1[i] = (char*) memalloc(MAX_SIZE);
      // for 2nd read, save extra fields for revComp(seq) and rev(qual)
      read2[i] = (char*) memalloc(MAX_SIZE);
    }
    char* header = (char*) memalloc(MAX_SIZE); // consensus header

    // process reads
    int len1 = 0, len2 = 0; // lengths of reads
    while (loadReads(in1, in2, read1, read2, header,
        &len1, &len2, offset, maxQual, gz1, gz2)) {

      // find optimal overlap
      float best = 1.0f;
      int pos = findPos(read1[SEQ], read2[SEQ + EXTRA + 1],
        read1[QUAL], read2[QUAL + EXTRA], len1, len2, overlap,
        dovetail, doveOverlap, mismatch, maxLen, &best);

      // print result
      if (pos == len1 - overlap + 1) {
        // stitch failure
        if (adaptOpt)
          printFail(out, out2, 1, log, 0, header, read1,
            read2, gzOut, lock + OUT, lock + LOG);
        else
          printFail(un1, un2, unOpt, log, logOpt, header,
            read1, read2, gzOut, lock + UN, lock + LOG);
      } else {
        // stitch success
        if (adaptOpt) {
          stitchRed += printResAdapt(out, out2, dove, doveOpt,
            header, read1, read2, len1, len2, pos, best,
            gzOut, lock);
        } else {
          printRes(out, log, logOpt, dove, doveOpt, aln, alnOpt,
            header, read1, read2, len1, len2, pos, best, offset,
            gzOut, fjoin, match, mism, lock);
          stitchRed++;
        }
      }

      count++;
    }

    // free memory
    free(header);
    for (int i = 0; i < FASTQ + EXTRA; i++) {
      if (i < FASTQ)
        free(read1[i]);
      free(read2[i]);
    }
    free(read1);
    free(read2);

  }  // END parallel

  // destroy omp locks
  for (int i = 0; i < 5; i++)
    omp_destroy_lock(&lock[i]);

  *stitch = stitchRed;
*/
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
int parseBAM(gzFile in, File out, bool gzOut, char* line,
    int chromLen, Chrom** chrom, int n_ref, int idx[],
    int* readLen, int* readMem, Read** unpaired, Read* dummy,
    unsigned long* totalLen, int* unmapped,
    int* paired, int* single, int* pairedPr,
    int* singlePr, int* supp, int* skipped, int* lowMapQ,
    int minMapQ, bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt) {

  // BAM fields to save
  int32_t refID, pos, l_seq, next_refID, next_pos, tlen;
  uint16_t n_cigar_op, flag;
  uint8_t mapq;
  uint32_t* cigar;
  uint8_t* seq;
  char* read_name, *qual, *extra;

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
    if (flag & 0xF00) {
      // skip supplementary/secondary alignments
      (*supp)++;
      continue;
    }
    Chrom* ref = chrom[idx[refID]];
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

    // save alignment information
    int length = calcDistBAM(l_seq, n_cigar_op, cigar); // distance to 3' end
    parseAlign(out, gzOut, readLen, readMem, &unpaired, dummy,
      read_name, flag, ref, pos,
      length, totalLen, paired, single, pairedPr, singlePr,
      singleOpt, extendOpt, extend, avgExtOpt);
    // NOTE: the following BAM fields are ignored:
    //   next_refID, next_pos, tlen, seq, qual, extra (optional fields)
  }

  // process single alignments w/ avgExtOpt
  if (avgExtOpt)
    *singlePr += saveAvgExt(*readLen, unpaired,
      *totalLen, *pairedPr);

  return count;
}

/* int readBAM()
 * Parse the header from a BAM file, then
 *   call parseBAM().
 */
int readBAM(gzFile in, File out, bool gzOut, char* line,
    unsigned long* totalLen,
    int* unmapped, int* paired, int* single,
    int* pairedPr, int* singlePr, int* supp, int* skipped,
    int* lowMapQ, int minMapQ, int xcount, char** xchrList,
    int* chromLen, Chrom*** chrom,
    int* readLen, int* readMem, Read** unpaired, Read* dummy,
    bool singleOpt, bool extendOpt, int extend, bool avgExtOpt) {
  // skip header
  int32_t l_text = readInt32(in, true);
  if (gzseek(in, l_text, SEEK_CUR) == -1)
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
      chromLen, chrom, xcount, xchrList);
  }

  return parseBAM(in, out, gzOut, line, *chromLen, *chrom,
    n_ref, idx, readLen, readMem, unpaired, dummy,
    totalLen, unmapped, paired, single, pairedPr,
    singlePr, supp, skipped, lowMapQ, minMapQ, singleOpt,
    extendOpt, extend, avgExtOpt);
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
void openFiles(char* outFile, File* out, File* out2,
    char* unFile, File* un1, File* un2,
    char* logFile, File* log,
    char* doveFile, File* dove, bool dovetail,
    char* alnFile, File* aln,
    bool gz) {

/*  if (adaptOpt) {
    if (interOpt)
      openWrite(outFile, out, gz);
    else if (! strcmp(outFile, "-"))
      exit(error("stdout + \"_1.fastq\"", ERROPENW));
    else if (! strcmp(outFile, "/dev/null")) {
      openWrite(outFile, out, gz);
      openWrite(outFile, out2, gz);
    } else {
      // add "_1.fastq" and "_2.fastq" extensions
      int add = strlen(ONEEXT) > strlen(TWOEXT) ?
        strlen(ONEEXT) + 1 : strlen(TWOEXT) + 1;
      char* outFile2 = memalloc(strlen(outFile) + add);
      strcpy(outFile2, outFile);
      strcat(outFile2, ONEEXT);
      openWrite(outFile2, out, gz);
      strcpy(outFile2, outFile);
      strcat(outFile2, TWOEXT);
      openWrite(outFile2, out2, gz);
      free(outFile2);
    }

  } else {*/
    openWrite(outFile, out, gz);

    // open optional files
/*    if (unFile != NULL) {
      if (interOpt)
        openWrite(unFile, un1, gz);
      else if (! strcmp(unFile, "-"))
        exit(error("stdout + \"_1.fastq\"", ERROPENW));
      else {
        // add "_1.fastq" and "_2.fastq" extensions
        int add = strlen(ONEEXT) > strlen(TWOEXT) ?
          strlen(ONEEXT) + 1 : strlen(TWOEXT) + 1;
        char* unFile2 = memalloc(strlen(unFile) + add);
        strcpy(unFile2, unFile);
        strcat(unFile2, ONEEXT);
        openWrite(unFile2, un1, gz);
        strcpy(unFile2, unFile);
        strcat(unFile2, TWOEXT);
        openWrite(unFile2, un2, gz);
        free(unFile2);
      }
    }*/
    if (logFile != NULL) {
      openWrite(logFile, log, false);
      fprintf(log->f, "Read\tOverlapLen\tStitchedLen\tMismatch\n");
    }
    if (alnFile != NULL)
      openWrite(alnFile, aln, false);
//  }

  if (dovetail && doveFile != NULL) {
    openWrite(doveFile, dove, false);
    fprintf(dove->f, "Read\tAdapter_R1\tAdapter_R2\n");
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
      exit(error("", ERRCLOSE));
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
    int skipped, int xcount, char** xchrList, int minMapQ,
    int lowMapQ, int paired, int orphan, int single,
    int singlePr, int pairedPr, unsigned long totalLen,
    bool singleOpt, bool extendOpt, int extend,
    bool avgExtOpt) {
  double avgLen = (double) totalLen / pairedPr;
  fprintf(stderr, "  SAM records analyzed: %10d\n", count);
  if (unmapped)
    fprintf(stderr, "    Unmapped:           %10d\n", unmapped);
  if (supp)
    fprintf(stderr, "    Supplementary:      %10d\n", supp);
  if (skipped) {
    fprintf(stderr, "    To skipped refs:    %10d\n", skipped);
    fprintf(stderr, "      (%s", xchrList[0]);
    for (int i = 1; i < xcount; i++)
      fprintf(stderr, ",%s", xchrList[i]);
    fprintf(stderr, ")\n");
  }
  if (lowMapQ)
    fprintf(stderr, "    MAPQ < %-2d:          %10d\n", minMapQ, lowMapQ);
  fprintf(stderr, "    Paired alignments:  %10d\n", paired);
  if (orphan)
    fprintf(stderr, "      orphan alignments:%10d\n", orphan);
  fprintf(stderr, "    Unpaired alignments:%10d\n", single);
  fprintf(stderr, "  BED intervals written:%10d\n", singlePr + pairedPr);
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

/* int checkOrphans()
 * Check for orphan paired alignments (missing mates).
 *   Free memory. Return count.
 */
int checkOrphans(Read* dummy) {
  int count = 0;
  Read* r = dummy->next;
  Read* tmp;
  while (r != NULL) {
    fprintf(stderr, "Warning! Read %s missing its pair -- not printed\n",
      r->name);
    count++;
    tmp = r->next;
    free(r->name);
    free(r);
    r = tmp;
  }
  dummy->next = NULL;
  return count;
}

/* void runProgram()
 * Controls the opening/closing of files,
 *   and analysis by readSAM() or readBAM().
 */
void runProgram(char* outFile, char* inFile, char* ctrlFile,
    bool singleOpt, bool extendOpt, int extend, bool avgExtOpt,
    int minMapQ,
    bool inter, char* unFile,
    char* logFile, int overlap, bool dovetail,
    char* doveFile, int doveOverlap, char* alnFile,
    int alnOpt, bool gzOut, bool fjoin,
    float mismatch, bool maxLen,
    int xcount, char** xchrList,
    bool verbose, int threads) {

  // open output files
  File out, out2, un1, un2, log, dove, aln; // output files
  openFiles(outFile, &out, &out2,
    unFile, &un1, &un2, logFile, &log,
    doveFile, &dove, dovetail, alnFile, &aln,
    gzOut);

  // initialize variables
  char* line = (char*) memalloc(MAX_SIZE);
  int chromLen = 0;       // number of reference sequences
  Chrom** chrom = NULL;   // array of reference sequences
  Read* dummy = (Read*) memalloc(sizeof(Read)); // head of linked list
  dummy->next = NULL;                           //   for paired alns
  unsigned long fragLen = 0;   // total length of treatment fragments

  // loop through input files (treatment and control)
  for (int i = 0; i < 2; i++) {
    char* file = (i ? ctrlFile : inFile);
    if (file == NULL) {
      savePileupNoCtrl(chrom, chromLen, fragLen);
      break;
    }

    char* end;
    char* filename = strtok_r(file, COM, &end);
    while (filename) {
      // open input files
      File in;
      bool gz = openRead(filename, &in);
      bool bam = checkBAM(in, gz);

      // process files
      if (verbose)
        fprintf(stderr, "Processing file: %s\n", filename);
      unsigned long totalLen = 0;  // total length of paired reads
      int unmapped = 0, paired = 0, single = 0, orphan = 0,
        pairedPr = 0, singlePr = 0, supp = 0, skipped = 0,
        lowMapQ = 0;  // counting variables
      Read** unpaired = NULL;       // for unpaired alns with
      int readLen = 0, readMem = 0; //   average-extension option
      int count;
      if (bam)
        count = readBAM(in.gzf, out, gzOut, line,
          &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ, minMapQ,
          xcount, xchrList, &chromLen, &chrom,
          &readLen, &readMem, unpaired, dummy,
          singleOpt, extendOpt, extend, avgExtOpt);
      else
        count = readSAM(in, gz, out, gzOut, line,
          un1, un2, unFile != NULL,
          log, logFile != NULL,
          overlap, dovetail, doveOverlap, dove,
          dovetail && doveFile != NULL, aln, alnOpt,
          mismatch, maxLen,
          &totalLen, &unmapped, &paired, &single,
          &pairedPr, &singlePr, &supp, &skipped, &lowMapQ, minMapQ,
          xcount, xchrList, &chromLen, &chrom,
          &readLen, &readMem, unpaired, dummy,
          singleOpt, extendOpt, extend, avgExtOpt,
          threads);

      // log counts
      orphan = checkOrphans(dummy);  // orphan paired alignments
      if (verbose)
        logCounts(count, unmapped, supp, skipped, xcount, xchrList,
          minMapQ, lowMapQ, paired, orphan, single, singlePr, pairedPr,
          totalLen, singleOpt, extendOpt, extend, avgExtOpt);

      // close input files
      if ( (gz && gzclose(in.gzf) != Z_OK) || (! gz && fclose(in.f)) )
        exit(error("", ERRCLOSE));

      filename = strtok_r(NULL, COM, &end);
    }

    // save pileup values
    fragLen = savePileup(chrom, chromLen, fragLen, i);

    // reset 'diff' array for each Chrom
    if (! i)
      for (int j = 0; j < chromLen; j++)
        if (chrom[j]->diff != NULL)
          for (int k = 0; k < chrom[j]->len; k++)
            chrom[j]->diff[k] = 0.0;

  }

  // save pileup values
  printPileup(out, gzOut, chrom, chromLen);

  // free memory
  if (xcount) {
    for (int i = 0; i < xcount; i++)
      free(xchrList[i]);
    free(xchrList);
  }
  for (int i = 0; i < chromLen; i++) {
    Chrom* chr = chrom[i];
    if (chr->ctrl != NULL) {
      free(chr->ctrl->end);
      free(chr->ctrl->cov);
      free(chr->ctrl);
    }
    if (chr->treat != NULL) {
      free(chr->treat->end);
      free(chr->treat->cov);
      free(chr->treat);
    }
    free(chr->diff);
    free(chr->name);
    free(chr);
  }
  free(chrom);
  free(dummy);
  free(line);

  // close files
  if ( ( gzOut && ( gzclose(out.gzf) != Z_OK ||
      (unFile != NULL && gzclose(un1.gzf) != Z_OK ) ) ) ||
      ( ! gzOut && ( fclose(out.f) ||
      (unFile != NULL && fclose(un1.f) ) ) ) ||
      (logFile != NULL && fclose(log.f)) ||
      (dovetail && doveFile != NULL && fclose(dove.f)) ||
      (alnFile != NULL && fclose(aln.f)) )
    exit(error("", ERRCLOSE));
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
  char* outFile = NULL, *inFile = NULL, *ctrlFile = NULL,
    *unFile = NULL, *logFile = NULL, *doveFile = NULL,
    *alnFile = NULL;
  char* xchrom = NULL;
  int extend = 0, minMapQ = 0,
    overlap = DEFOVER, doveOverlap = DEFDOVE,
    threads = DEFTHR;
  float mismatch = DEFMISM;
  bool singleOpt = false, extendOpt = false, avgExtOpt = false;
  bool gzOut = false;
  bool dovetail = false, maxLen = true,
    fjoin = false,
    verbose = false;

  // parse argv
  int c;
  while ( (c = getopt_long(argc, argv, OPTIONS, long_options, NULL)) != -1 )
    switch (c) {
      case INFILE: inFile = optarg; break;
      case CTRLFILE: ctrlFile = optarg; break;
      case OUTFILE: outFile = optarg; break;
      case SINGLEOPT: singleOpt = true; break;
      case EXTENDOPT: extend = getInt(optarg); extendOpt = true; break;
      case AVGEXTOPT: avgExtOpt = true; break;
      case XCHROM: xchrom = optarg; break;
      case MINMAPQ: minMapQ = getInt(optarg); break;
      case GZOPT: gzOut = true; break;

      case VERBOSE: verbose = true; break;
      case VERSOPT: printVersion(); break;
      case HELP: usage(); break;

      case DOVEOPT: dovetail = true; break;
      case FJOINOPT: fjoin = true; break;
      case UNFILE: unFile = optarg; break;
      case LOGFILE: logFile = optarg; break;
      case DOVEFILE: doveFile = optarg; break;
      case ALNFILE: alnFile = optarg; break;
      case OVERLAP: overlap = getInt(optarg); break;
      case DOVEOVER: doveOverlap = getInt(optarg); break;
      case MISMATCH: mismatch = getFloat(optarg); break;
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

  // save list of chromosomes to ignore
  int xcount = 0;
  char** xchrList = NULL;
  if (xchrom != NULL)
    xcount = saveXChrom(xchrom, &xchrList);

  bool inter = false;  // interleaved input
  if (overlap <= 0 || doveOverlap <= 0)
    exit(error("", ERROVER));
  if (mismatch < 0.0f || mismatch >= 1.0f)
    exit(error("", ERRMISM));
  if (threads < 1)
    exit(error("", ERRTHREAD));

  // adjust parameters for adapter-removal mode
  if (true) {
    dovetail = true;
    unFile = logFile = alnFile = NULL;
  }
  int alnOpt = (alnFile != NULL ? 1 : 0);

  // send arguments to runProgram()
  runProgram(outFile, inFile, ctrlFile,
    singleOpt, extendOpt, extend, avgExtOpt, minMapQ,
    inter, unFile,
    logFile, overlap, dovetail, doveFile, doveOverlap,
    alnFile, alnOpt, gzOut, fjoin,
    mismatch, maxLen,
    xcount, xchrList,
    verbose, threads);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getArgs(argc, argv);
  return 0;
}
