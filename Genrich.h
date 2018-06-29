/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.0
*/
#define VERSION     "0.0"

// constants
#define MAX_SIZE    1024    // maximum length of input lines (incl. seq/qual)
#define NOTMATCH    1.5f    // stitch failure
#define COM         ", "    // separator for input file names
#define CSV         ",\t"   // separator for quality score profile
#define NA          "NA"    // n/a (for output log file)
#define TAB         "\t"

// default parameter values
#define DEFOVER     20      // min. overlap
#define DEFDOVE     50      // min. overlap for dovetailed alignments
#define DEFMISM     0.1f    // mismatch fraction
#define OFFSET      33      // fastq quality offset (Sanger = 33)
#define MAXQUAL     40      // maximum quality score (0-based)
#define DEFTHR      1       // number of threads

// fastq parts
enum fastq { HEAD, SEQ2, PLUS, QUAL2, FASTQ };  // lines of a fastq read
#define BEGIN       '@'     // beginning of header line
#define PLUSCHAR    '+'     // beginning of 3rd line
#define EXTRA       2       // save 2 extra strings for 2nd read:
                            //   revComp(seq) and rev(qual)

// SAM fields
enum sam { NIL, QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
  PNEXT, TLEN, SEQ, QUAL };


// command-line options
#define OPTIONS     "hi:o:ya:xc:f:l:m:p:de:C:sj:bzgq:n:vV"
#define HELP        'h'
#define INFILE      'i'
#define OUTFILE     'o'
#define SINGLEOPT   'y'
#define EXTENDOPT   'a'
#define AVGEXTOPT   'x'
#define XCHROM      'c'
#define UNFILE      'f'
#define LOGFILE     'l'
#define OVERLAP     'm'
#define MISMATCH    'p'
#define DOVEOPT     'd'
#define DOVEOVER    'e'
#define DOVEFILE    'C'
#define MAXOPT      's'
#define ALNFILE     'j'
#define DIFFOPT     'b'
#define GZOPT       'z'
#define FJOINOPT    'g'
#define QUALITY     'q'
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
#define ONEEXT      "_1.fastq"
#define TWOEXT      "_2.fastq"
#define GZEXT       ".gz"   // for gzip compression

// OMP locks
enum omp_locks { OUT, UN, LOG, DOVE, ALN, OMP_LOCKS };

// error messages
enum errCode { ERRFILE, ERROPEN, ERRCLOSE, ERROPENW, ERRUNK,
  ERRMEM, ERRSEQ, ERRQUAL, ERRHEAD, ERRINT, ERRFLOAT, ERRPARAM,
  ERROVER, ERRMISM, ERRINFO, ERRSAM, ERRREP, ERRCHROM, ERREXTEND,
  ERROFFSET, ERRUNGET, ERRGZIP,
  ERRTHREAD, ERRNAME, ERRRANGE, ERRDEFQ, ERRCIGAR, DEFERR
};
const char* errMsg[] = { "Need input/output files",
  ": cannot open file for reading",
  "Cannot close file",
  ": cannot open file for writing",
  ": unknown nucleotide",
  "Cannot allocate memory",
  "Cannot load sequence",
  "Sequence/quality scores do not match",
  ": not matched in input files",
  ": cannot convert to int",
  ": cannot convert to float",
  ": unknown command-line argument",
  "Overlap must be greater than 0",
  ": mismatch between sequence length and CIGAR",
  ": no sequence information (SEQ or CIGAR)",
  ": poorly formatted SAM record",
  ": read has repeated information in SAM",
  ": cannot find reference sequence name in SAM header",
  "Extension length must be >= 0",

  ": quality score outside of set range",
  "Failure in ungetc() call",
  "Cannot pipe in gzip compressed file (use zcat instead)",
  "Number of threads must be >= 1",
  ": output filename cannot start with '-'",
  ": file missing values for quality score range",
  "Cannot increase max. quality score with default error profile",
  ": unknown Op in CIGAR",
  "Unknown error"
};

// generic File type
typedef union file {
  FILE* f;
  gzFile gzf;
} File;

typedef struct chrom {
  char* name;
  int len;
  struct chrom* next;
} Chrom;

typedef struct read {
  int posR1;
  int posR2;
  bool strand;
  char* name;
  char* chrom;
  struct read* next;
//  struct read* prev;
} Read;
