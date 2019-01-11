/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  June 2018

  Finding sites of enrichment from genome-wide assays.

  Version 0.5
*/
#define VERSION     "0.5"

// macros
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// constants
#define MAX_SIZE    65520   // maximum length of input SAM/BAM alignments
#define MAX_ALNS    128     // maximum number of alignments per read/pair
                            //   - also used as max. read name length,
                            //     and for various dynamic memory allocs
#define HASH_SIZE   1310417 // size of hashtable for p-values
#define TAB         "\t"    // separator for SAM/BED fields
#define TABN        "\t\n"  // separator for final BED field
#define COL         ":"     // separator for SAM optional fields (TAG:TYPE:VALUE)
#define COM         ", "    // separator for input file names / ref. names
#define NA          "NA"    // results not available
#define GZEXT       ".gz"   // extension for gzip-compressed files
#define SKIP        -1.0f   // stats for a genomic region to be skipped

// default parameter values
#define DEFQVAL     0.05f   // default q-value
#define DEFAUC      20.0f   // area under the curve for peak calling
#define DEFMAXGAP   100     // maximum gap between significant sites
#define DEFMINLEN   0       // minimum length of a peak
#define DEFATAC     100     // interval length for ATAC-seq mode

// SAM fields
enum sam { NIL, QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
  PNEXT, TLEN, SEQ, QUAL };
#define SAMQUAL     33        // quality score offset
#define SCORE       "AS"      // extra field: alignment score
#define NOSCORE     -FLT_MAX  // for alignments with no alignment score(s)

// alignment types
enum alignType { PAIRED, SINGLE, DISCORD };

// fields of a bedGraph file
enum bedGraph { CHR, START, END };

// constants for log-normal p-value calculation
#define LOGSQRT     0.445999019652555   // log(sqrt(2.44))
#define SQRTLOG     0.944456478248262   // sqrt(log(2.44))

// command-line options
#define OPTIONS     "ht:c:o:f:k:b:zyw:xjd:e:E:m:s:p:q:a:l:g:rR:XPvV"
#define HELP        'h'
#define INFILE      't'
#define CTRLFILE    'c'
#define OUTFILE     'o'
#define LOGFILE     'f'
#define PILEFILE    'k'
#define BEDFILE     'b'
#define GZOPT       'z'
#define UNPAIROPT   'y'
#define EXTENDOPT   'w'
#define AVGEXTOPT   'x'
#define ATACOPT     'j'
#define ATACLEN     'd'
#define XCHROM      'e'
#define XFILE       'E'
#define MINMAPQ     'm'
#define ASDIFF      's'
#define PVALUE      'p'
#define QVALUE      'q'
#define MINAUC      'a'
#define MINLEN      'l'
#define MAXGAP      'g'
#define DUPSOPT     'r'
#define DUPSFILE    'R'
#define NOPEAKS     'X'
#define PEAKSONLY   'P'
#define VERBOSE     'v'
#define VERSOPT     'V'

static struct option long_options[] = {
  {"help", no_argument, NULL, HELP},
  {"verbose", no_argument, NULL, VERBOSE},
  {"version", no_argument, NULL, VERSOPT},
  {0, 0, 0, 0}
};

// error messages
enum errCode { ERRFILE, ERROPEN, ERROPENW, ERRCLOSE,
  ERRMEM, ERRINT, ERRFLOAT, ERRPARAM, ERREXTEND, ERRATAC,
  ERRPQVAL, ERRASDIFF, ERRMINAUC, ERRMINLEN, ERRMISM,
  ERRINFO, ERRSAM, ERRCHROM, ERRHEAD, ERRBAM, ERRGEN,
  ERREXPT, ERRCHRLEN, ERRCTRL, ERRPOS, ERRSORT, ERRTYPE,
  ERRAUX, ERRBED, ERRLINEAR, ERRINDEX, ERRLOGIDX, ERRLOG,
  ERRISSUE, ERRALNS, ERRPILE, ERRPVAL, ERRARR, ERRARRC,
  ERRDF, ERRALNTYPE, ERRUNGET, ERRGZIP, ERRNAME, ERRCIGAR,
  DEFERR
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
  "Minimum AUC must be >= 0.0",
  "Minimum peak length must be >= 0",
  ": mismatch between sequence length and CIGAR",
  ": no sequence information (SEQ or CIGAR)",
  ": poorly formatted SAM/BAM record",
  ": cannot find reference sequence name in SAM header",
  ": misplaced SAM header line",
  "Cannot parse BAM file",
  "No analyzable genome (length=0)",
  "Experimental sample has no analyzable fragments",
  ": reference sequence has different lengths in BAM/SAM files",
  ": reference sequence missing from control sample(s)",
  ": read aligned beyond reference end",
  "SAM/BAM file not sorted by queryname (samtools sort -n)",
  ": unknown value type in BAM auxiliary field",
  "Poorly formatted BAM auxiliary field",
  ": poorly formatted BED record",
  "Linear template with >2 reads -- not allowed",
  "Unknown index of paired alignment",
  ": cannot find field in header of bedgraph-ish log file",
  "Poorly formatted bedgraph-ish log record",
  "\n  (internal error: please open an Issue on https://github.com/jsh58/Genrich)",
  "Disallowed number of alignments",
  "Invalid pileup value (< 0)",
  "Failure collecting p-values",
  "Failure creating experimental pileup",
  "Failure creating control pileup",
  "Invalid df in pchisq()",
  "Invalid alignment type",
  "Failure in ungetc() call",
  "Cannot pipe in gzip-compressed file (use zcat instead)",
  ": output filename cannot start with '-'",
  ": unknown Op in CIGAR",
  "Unknown error"
};

// generic File type
typedef union file {
  FILE* f;
  gzFile gzf;
} File;

typedef struct bed {
  uint32_t pos[2];
  char* name;     // chromosome name
} Bed;

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
  bool save;          // chromosome to be saved? (by sample)
  uint32_t* bed;      // coordinates (paired) of regions to be excluded
  int bedLen;         // number of coordinates of regions to be excluded
  Diff* diff;         // arrays for keeping track of pileup changes
  Pileup* expt;       // pileup arrays for experimental sample(s)
  uint32_t exptLen;   // length of pileup arrays for experimental sample(s) (dynamic)
  uint32_t exptMem;   // length of pileup arrays for experimental sample(s) (in memory)
  Pileup* ctrl;       // pileup arrays for control sample(s)
  uint32_t ctrlLen;   // length of pileup arrays for control sample(s) (dynamic)
  uint32_t ctrlMem;   // length of pileup arrays for control sample(s) (in memory)
  Pileup** pval;      // "pileup" arrays for p-values
  uint32_t* pvalLen;  // lengths of "pileup" arrays for p-values
  uint8_t sample;     // count of samples with p-value arrays saved
  Pileup* qval;       // "pileup" arrays for q-values
} Chrom;

typedef struct aln {
  uint32_t pos[2];  // positions of the alignment
  float score;      // alignment score (sum of scores for paired alns)
  bool primary;     // primary alignment?
  bool paired;      // properly paired alignment?
  bool full;        // both parts of paired aln analyzed? (only for paired alns)
  bool first;       // which read of a pair this is (true -> R1; false -> R2)
  bool strand;      // which strand aln is on (only for unpaired alns)
  uint8_t count;    // value of aln (only for unpaired alns with avg-ext option)
  char* name;       // read name (only for unpaired alns with avg-ext option)
  Chrom* chrom;     // reference sequence
} Aln;

typedef struct hashAln {
  char* name;       // read name
  Chrom* chrom;     // reference sequence
  Chrom* chrom1;    // other reference sequence (discordant only)
  uint32_t pos;     // position of the alignment
  uint32_t pos1;    // other position of the aln (paired and discordant only)
  bool strand;      // strand of aln (discordant and singletons only)
  bool strand1;     // other strand of aln (discordant only)
  struct hashAln* next;
} HashAln;

typedef struct read {
  char* name;       // read name
  Aln* aln;         // array of alignments
  uint8_t alnLen;   // length of alignment array
  Aln* alnR2;       // array of alignments for R2 (discordant alns only)
  uint8_t alnLenR2; // length of alnR2 array (discordant alns only)
  uint16_t qual;    // sum of quality scores
  bool first;       // true -> R1; false -> R2 (singleton alns only)
  float score;      // min. alignment score
  float scoreR2;    // min. alignment score for R2 (discordant alns only)
} Read;
