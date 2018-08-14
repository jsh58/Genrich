/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.1
*/
#define VERSION     "0.1"

// macros
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// constants
#define MAX_SIZE    65520   // maximum length of input SAM/BAM alignments
#define MAX_ALNS    100     // maximum number of alignments per read/pair
#define TAB         "\t"    // separator for SAM fields
#define COM         ", "    // separator for input file names / ref. names

// default parameter values
#define DEFQVAL     0.05    // default q-value
#define DEFMINLEN   100     // minimum length of a peak
#define DEFMAXGAP   100     // maximum gap between significant sites
#define DEFTHR      1       // number of threads

// SAM fields
enum sam { NIL, QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
  PNEXT, TLEN, SEQ, QUAL };

// command-line options
#define OPTIONS     "ht:c:o:b:zya:xe:m:p:q:g:l:n:vV"
#define HELP        'h'
#define INFILE      't'
#define CTRLFILE    'c'
#define OUTFILE     'o'
#define LOGFILE     'b'
#define GZOPT       'z'
#define SINGLEOPT   'y'
#define EXTENDOPT   'a'
#define AVGEXTOPT   'x'
#define XCHROM      'e'
#define MINMAPQ     'm'
#define PVALUE      'p'
#define QVALUE      'q'
#define MAXGAP      'g'
#define MINLEN      'l'

#define THREADS     'n'
#define VERBOSE     'v'
#define VERSOPT     'V'

static struct option long_options[] = {
  {"help", no_argument, NULL, HELP},
  {"verbose", no_argument, NULL, VERBOSE},
  {"version", no_argument, NULL, VERSOPT},
  {0, 0, 0, 0}
};

// extensions for output files
#define GZEXT       ".gz"   // for gzip compression

// OMP locks
enum omp_locks { OUT, UN, LOG, DOVE, ALN, OMP_LOCKS };

// error messages
enum errCode { ERRFILE, ERROPEN, ERROPENW, ERRCLOSE, ERRMEM,
  ERRSEQ, ERRQUAL, ERRHEAD, ERRINT, ERRFLOAT, ERRPARAM,
  ERRMISM, ERRINFO, ERRSAM, ERRREP, ERRCHROM, ERREXTEND,
  ERRBAM, ERRGEN, ERRCHRLEN, ERRCTRL, ERRPOS, ERRSORT,
ERRUNGET, ERRGZIP,
  ERRTHREAD, ERRNAME, ERRCIGAR, DEFERR
};
const char* errMsg[] = { "Need input/output files",
  ": cannot open file for reading",
  ": cannot open file for writing",
  ": cannot close file",
  "Cannot allocate memory",
  "Cannot load sequence",
  "Sequence/quality scores do not match",
  ": not matched in input files",
  ": cannot convert to int",
  ": cannot convert to float",
  ": unknown command-line argument",
  ": mismatch between sequence length and CIGAR",
  ": no sequence information (SEQ or CIGAR)",
  ": poorly formatted SAM/BAM record",
  ": read has repeated information in SAM",
  ": cannot find reference sequence name in SAM header",
  "Extension length must be >= 0",
  "Cannot parse BAM file",
  "No analyzable genome (length=0)",
  ": reference sequence has different lengths in BAM/SAM files",
  ": reference sequence missing from control sample(s)",
  ": read aligned beyond reference end",
  "SAM/BAM file not sorted by queryname (samtools sort -n)",

  "Failure in ungetc() call",
  "Cannot pipe in gzip compressed file (use zcat instead)",
  "Number of threads must be >= 1",
  ": output filename cannot start with '-'",
  ": unknown Op in CIGAR",
  "Unknown error"
};

// generic File type
typedef union file {
  FILE* f;
  gzFile gzf;
} File;

typedef struct pileup {
  unsigned int* end;
  float* cov;
} Pileup;

typedef struct chrom {
  char* name;
  uint32_t len;
  bool skip;
  float* diff;
  Pileup* treat;
  int treatLen; // length of pileup arrays for treatment sample(s)
  Pileup* ctrl;
  int ctrlLen;  // length of pileup arrays for control sample(s)
  Pileup* pval;
  Pileup* qval;
  int pvalLen;  // length of pileup arrays for p- and q-values
} Chrom;

typedef struct aln {
  uint32_t pos[2];
  bool strand;  // only for SE alignments
  bool paired;  // PE alignment
  bool full;    // both parts of PE aln analyzed
  float val;    // value for aln (only for singletons with avg-ext option)
  char* name;   // read name (only for singletons with avg-ext option)
  Chrom* chrom;
} Aln;

typedef struct read {
  char* name;
  Aln** aln;
  int alnLen;
  struct read* next;
} Read;
