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
#define MAX_ALNS    128     // maximum number of alignments per read/pair
                            //   - also used as max. read name length,
                            //     and for various dynamic memory allocs
#define HASH_SIZE   1310417 // size of hashtable for p-values
#define TAB         "\t"    // separator for SAM fields
#define COL         ":"     // separator for SAM optional fields (TAG:TYPE:VALUE)
#define COM         ", "    // separator for input file names / ref. names

// default parameter values
#define DEFQVAL     0.05    // default q-value
#define DEFMINLEN   100     // minimum length of a peak
#define DEFMAXGAP   100     // maximum gap between significant sites
#define DEFATAC     100     // interval length for ATAC-seq mode
#define DEFTHR      1       // number of threads

// SAM fields
enum sam { NIL, QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
  PNEXT, TLEN, SEQ, QUAL };
#define SCORE       "AS"      // extra field: alignment score
#define NOSCORE     -FLT_MAX  // for alignments with no alignment score(s)

// command-line options
#define OPTIONS     "ht:c:o:b:zya:xjd:e:m:s:p:q:g:l:n:vV"
#define HELP        'h'
#define INFILE      't'
#define CTRLFILE    'c'
#define OUTFILE     'o'
#define LOGFILE     'b'
#define GZOPT       'z'
#define SINGLEOPT   'y'
#define EXTENDOPT   'a'
#define AVGEXTOPT   'x'
#define ATACOPT     'j'
#define ATACLEN     'd'
#define XCHROM      'e'
#define MINMAPQ     'm'
#define ASDIFF      's'
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
enum errCode { ERRFILE, ERROPEN, ERROPENW, ERRCLOSE,
  ERRMEM, ERRINT, ERRFLOAT, ERRPARAM, ERREXTEND, ERRATAC,
  ERRPQVAL, ERRASDIFF, ERRMISM, ERRINFO, ERRSAM, ERRCHROM,
  ERRHEAD, ERRBAM, ERRGEN, ERRTREAT, ERRCHRLEN, ERRCTRL,
  ERRPOS, ERRSORT, ERRTYPE, ERRAUX, ERRISSUE, ERRALNS,
  ERRPILE, ERRPVAL, ERRARR, ERRARRC, ERRUNGET, ERRGZIP,
  ERRTHREAD, ERRNAME, ERRCIGAR, DEFERR
};
const char* errMsg[] = { "Need input/output files",
  ": cannot open file for reading",
  ": cannot open file for writing",
  ": cannot close file",
  "Cannot allocate memory",
  ": cannot convert to int",
  ": cannot convert to float",
  ": unknown command-line argument",
  "Extension length must be > 0",
  "ATAC-seq interval length must be > 0",
  "p-/q-value must be in (0,1]",
  "Secondary alignment score threshold must be >= 0.0",

  ": mismatch between sequence length and CIGAR",
  ": no sequence information (SEQ or CIGAR)",
  ": poorly formatted SAM/BAM record",
  ": cannot find reference sequence name in SAM header",
  ": misplaced SAM header line",
  "Cannot parse BAM file",
  "No analyzable genome (length=0)",
  "Treatment sample(s) have no analyzable fragments",
  ": reference sequence has different lengths in BAM/SAM files",
  ": reference sequence missing from control sample(s)",
  ": read aligned beyond reference end",
  "SAM/BAM file not sorted by queryname (samtools sort -n)",
  ": unknown value type in BAM auxiliary field",
  "Poorly formatted BAM auxiliary field",
  "\n  (internal error: please open an Issue on https://github.com/jsh58/Genrich)",
  "Disallowed number of alignments",
  "Invalid pileup value (< 0)",
  "Failure collecting p-values",
  "Failure creating treatment pileup",
  "Failure creating control pileup",

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

typedef struct hash {
  float val;      // p-value
  uint64_t len;   // length of genome with that p-value
  struct hash* next;
} Hash;

typedef struct pileup {
  uint32_t* end;  // array of end coordinates
  float* cov;     // array of pileup values
} Pileup;

typedef struct diff {
  uint8_t* frac;  // fractions of a count (8-bit encoded)
  int16_t* cov;   // int counts
} Diff;

typedef struct chrom {
  char* name;         // name of chromosome (reference sequence)
  uint32_t len;       // length of chromosome
  bool skip;          // chromosome to be skipped?
  Diff* diff;         // arrays for keeping track of pileup changes
  Pileup* treat;      // pileup arrays for treatment sample(s)
  uint32_t treatLen;  // length of pileup arrays for treatment sample(s)
  Pileup* ctrl;       // pileup arrays for control sample(s)
  uint32_t ctrlLen;   // length of pileup arrays for control sample(s)
  Pileup* pval;       // "pileup" arrays for p-values
  Pileup* qval;       // "pileup" arrays for q-values
  uint32_t pvalLen;   // length of "pileup" arrays for p- and q-values
} Chrom;

typedef struct aln {
  uint32_t pos[2];  // positions of the alignment
  float score;      // alignment score (sum of scores for paired alns)
  bool primary;     // primary alignment?
  bool paired;      // properly paired alignment?
  bool full;        // both parts of paired aln analyzed? (only for paired alns)
  bool strand;      // which strand aln is on (only for singleton alignments)
  bool first;       // which read of a pair this is (only for singleton alignments)
  uint8_t count;    // value of aln (only for singletons with avg-ext option)
  char* name;       // read name (only for singletons with avg-ext option)
  Chrom* chrom;     // reference sequence
} Aln;

/*
typedef struct read {
  char* name;
  Aln** aln;
  int alnLen;
  struct read* next;
} Read;
*/
