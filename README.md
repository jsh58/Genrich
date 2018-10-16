# Genrich: detecting sites of genomic enrichment

## Table of Contents
* [Introduction](#intro)
  * [Quick start](#quick)
  * [Software compilation](#compile)
  * [Usage message](#usage)
* [Attributes](#attributes)
  * [Peak-calling method](#method)
  * [Alignment parsing](#alignment)
  * [Multiple replicates](#replicate)
  * [Multimapping reads](#multimap)
  * [Genome length calculation](#genomelen)
  * [Control/background pileup calculation](#pileup)
  * [*p*-value calculation](#pvalue)
  * [*q*-value calculation](#qvalue)
* [I/O files and options](#files)
  * [Required files](#required)
  * [Optional files](#optional)
* [Filtering options](#filter)
  * [Unpaired alignments](#unpaired)
* [ATAC-seq mode](#atacseq)
* [Peak-calling parameters](#peakcalling)
* [Miscellaneous](#misc)
  * [Warning messages](#warning)
* [Contact](#contact)
<br><br>

## Introduction <a name="intro"></a>

Genrich is a peak-caller for genomic enrichment assays (e.g. ChIP-seq, ATAC-seq).  It analyzes alignment files generated following the assay and produces a file detailing peaks of significant enrichment.


### Quick start <a name="quick"></a>

Given:
* `sample.bam` (alignment file)
* `Genrich` (compiled as described [below](#compile))

To produce a file listing regions of genomic enrichment:
```
$ ./Genrich  -t sample.bam  -o sample.narrowPeak  -v
```

### Software compilation <a name="compile"></a>

The software can be downloaded from [GitHub](https://github.com/jsh58/Genrich).

A Makefile is provided for compilation with [GCC](https://gcc.gnu.org/releases.html), and [zlib](http://zlib.net) is also required.  The program has been tested after compilation with GCC 5.4.0 and zlib 1.2.8.

To compile, run `make` in the folder in which the software was downloaded.  The executable `Genrich` should be produced.


### Usage message <a name="usage"></a>

```
Usage: ./Genrich  -t <file>  -o <file>  [optional arguments]
Required arguments:
  -t  <file>       Input SAM/BAM file(s) for treatment sample(s)
  -o  <file>       Output peak file (in ENCODE narrowPeak format)
Optional I/O arguments:
  -c  <file>       Input SAM/BAM file(s) for control sample(s)
                     (matched with -t files; 'null' if missing)
  -f  <file>       Output bedgraph-ish file for p/q values
  -k  <file>       Output bedgraph-ish file for pileups and p-values
  -b  <file>       Output BED file for reads/fragments/intervals
Filtering options:
  -e  <arg>        Comma-separated list of chromosomes to ignore
  -E  <file>       Input BED file(s) of genomic regions to ignore
  -m  <int>        Minimum MAPQ to keep an alignment (def. 0)
  -s  <float>      Keep sec alns with AS >= bestAS - <float> (def. 0)
  -y               Keep unpaired alignments (def. false)
  -w  <int>        Keep unpaired alns, lengths changed to <int>
  -x               Keep unpaired alns, lengths changed to paired avg
Options for ATAC-seq:
  -j               Use ATAC-seq mode (def. false)
  -d  <int>        Expand cut sites to <int> bp (def. 100)
Options for peak calling:
  -q  <float>      Maximum q-value (FDR-adjusted p-value; def. 0.05)
  -p  <float>      Maximum p-value (overrides -q if set)
  -a  <float>      Minimum AUC for a peak (def. 20.0)
  -l  <int>        Minimum length of a peak (overrides -a if set)
  -g  <int>        Maximum distance between signif. sites (def. 100)
Other options:
  -z               Option to gzip-compress output(s)
  -v               Option to print status updates/counts to stderr
```

## Attributes <a name="attributes"></a>

### Peak-calling method <a name="method"></a>

Here is an overview of the method used by Genrich to identify peaks (Fig. 1):
* Parse alignments for the treatment sample and create a treatment "pileup" by counting the DNA fragments that cover each position of the genome (additional information about alignment parsing can be found [here](#alignment)).
* Create a control pileup using the control sample (if available) and background level (additional information about control/background pileup calculation can be found [here](#pileup)).
* Calculate *p*-values for each genomic position, as described [here](#pvalue).
* Convert *p*-values to *q*-values, as described [here](#qvalue).
* Calculate the "area under the curve" (AUC) for all regions whose -log(*q*) values rise above the statistical threshold.
* Combine nearby regions and call peaks whose total AUC is above a threshold (details of peak-calling parameters can be found [here](#peakcalling)).

<figure>
  <img src="figures/figure1.png" alt="Peak-calling by Genrich" width="800">
  <figcaption><strong>Figure 1.  Peak-calling by Genrich.</strong>  Visualization by <a href="http://software.broadinstitute.org/software/igv/">IGV</a>.</figcaption>
</figure>
<br><br>

### Alignment parsing <a name="alignment"></a>

Genrich analyzes paired-end reads aligned to a reference genome.  It correctly infers full fragments as spanning between the 5' ends of two properly paired alignments.  By default, it does not consider unpaired ("singleton") alignments, although there are three options for keeping such alignments, as described [here](#unpaired).

An alternative analysis mode for ATAC-seq is also provided by Genrich, as described [here](#atacseq).


### Multiple replicates <a name="replicate"></a>

Genrich calls peaks for multiple replicates collectively.  First, it analyzes the replicates separately, with [*p*-values calculated](#pvalue) for each. At each genomic position, the multiple replicates' *p*-values are then combined by [Fisher's method](https://en.wikipedia.org/wiki/Fisher's_method#Application_to_independent_test_statistics).  The combined *p*-values are [converted to *q*-values](#qvalue), and peaks are called.  This obviates the need for [IDR](https://www.encodeproject.org/software/idr/) (you're welcome!).


### Multimapping reads <a name="multimap"></a>

Genrich analyzes reads/fragments that map to multiple locations in the genome by adding a fractional count to each location.  This allows for peak detection in regions of the genome that are otherwise inaccessible to the assay.
* The input SAM/BAM file(s) must list secondary alignments for multimapping reads/fragments.
  * The short read aligner [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) produces secondary alignments in either [`-k <int>` mode](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#k-mode-search-for-one-or-more-alignments-report-each) or [`-a` mode](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#a-mode-search-for-and-report-all-alignments).
  * The short read aligner [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) does not produce secondary alignments.
* To avoid excessive memory usage and the imprecision inherent in floating-point values, a maximum of 10 alignments per read is analyzed by Genrich.
* Additional information can be found in the description of the [`-s` parameter](#sparam).


### Genome length calculation <a name="genomelen"></a>

Genrich computes the genome length as the sum of the lengths of the chromosomes (reference sequences) in the header of the SAM/BAM file.  The length is reduced if the user specifies chromosomes (`-e`) or genomic regions (`-E`) to be ignored, as described [below](#filter).  The calculated length is used for calculating a [background pileup value](#pileup) and [*q*-values](#qvalue).


### Control/background pileup calculation <a name="pileup"></a>

The background pileup value is calculated by dividing the total sequence information in the treatment sample by the [calculated genome length](#genomelen).  The net control pileup value at a particular genomic position is the maximum of the background pileup value and the pileup of the control sample at that position (if a control sample is specified).  Note that control pileups are scaled to match the treatment, based on the total sequence information in each (the scaling factor is given as part of the verbose [`-v`] output).


### *p*-value calculation <a name="pvalue"></a>

The *p*-values are calculated for each base of the genome assuming an [exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution#Alternative_parameterization) with the control/background pileup value as the parameter *&beta;*.  The exponential distribution has two attributes that make it suitable for this calculation:
* Because it is a continuous probability distribution, fractional treatment pileups can be considered.  This is important for reads/fragments with [multiple alignments](#multimap).
* The *p*-value for a treatment value *x* and parameter *&beta;* is <i>e<sup> -x / &beta;</sup></i>.  Therefore, multiplying *x* and *&beta;* by the same factor results in the same *p*-value, which implies that the *p*-values are robust to sequencing depth artifacts.


### *q*-value calculation <a name="qvalue"></a>

The *q*-value for each base of the genome is calculated from the *p*-value using the [Benjamini-Hochberg procedure](http://www.math.tau.ac.il/~ybenja/MyPapers/benjamini_hochberg1995.pdf).  Note that the [calculated genome length](#genomelen) is used as the number of hypothesis tests (*m*).


## I/O files and options <a name="files"></a>

### Required files <a name="required"></a>

```
  -t  <file>       Input SAM/BAM file(s) for treatment sample(s)
```
* Genrich analyzes alignment files in [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).  SAM files must have a header.
* SAM/BAM files for [multiple replicates](#replicate) can be specified, comma-separated (or space-separated, in quotes).
* Multiple SAM/BAM files for a single replicate should be combined in advance via `samtools merge`.
* The SAM/BAM files should be name sorted (via `samtools sort -n`).  As of [Version 0.3](https://github.com/jsh58/Genrich/releases/tag/v0.3), unsorted SAM/BAM files are allowed, but this is likely to change.
* Genrich will read from `stdin` with `-t -`.
<br><br>

```
  -o  <file>       Output peak file (in ENCODE narrowPeak format)
```
* As indicated, the output file is in [ENCODE narrowPeak format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12).  Here are additional details of the fields:
<table>
  <tr>
    <td align="center">4. name</td>
    <td><code>peak_N</code>, where <code>N</code> is the 0-based count</td>
  </tr>
  <tr>
    <td align="center">5. score</td>
    <td>10*qValue (or 10*pValue if <code>-p</code> is set), rounded to the nearest int (max. 1000)</td>
  </tr>
  <tr>
    <td nowrap align="center">7. signalValue</td>
    <td>Total area under the curve, in default peak-calling mode.  In minimum-length peak-calling mode (<code>-l</code>), the summit fold-enrichment (treatment / control) with one replicate, or the qValue (pValue with <code>-p</code>) with multiple replicates.</td>
  </tr>
  <tr>
    <td align="center">8. pValue</td>
    <td>Summit -log<sub>10</sub>(<i>p</i>-value)</td>
  </tr>
  <tr>
    <td align="center">9. qValue</td>
    <td>Summit -log<sub>10</sub>(<i>q</i>-value), or -1 with <code>-p</code></td>
  </tr>
  <tr>
    <td align="center">10. peak</td>
    <td>Summit position: the midpoint of the peak interval reaching the highest significance (the longest interval in case of ties)</td>
  </tr>
</table>

* Here is the part of the output file corresponding to the peaks called in Figure 1:
```
chr1    1565272    1565335    peak_253     43    .      92.013634     6.853281     4.313010     25
chr1    1565608    1566028    peak_254    129    .    1473.275024    15.990990    12.873972    259
```
* Genrich will write to `stdout` with `-o -`.


### Optional files <a name="optional"></a>

```
  -c  <file>       Input SAM/BAM file(s) for control sample(s)
                     (matched with -t files; 'null' if missing)
```
* Alignment files for control samples can be specified.  As indicated, they should be matched with treatment files.
* SAM/BAM files for multiple replicates can be specified, comma-separated (or space-separated, in quotes).  Missing control files should be indicated with `null`.
<br><br>

```
  -f  <file>       Output bedgraph-ish file for p/q values
```
* With a single replicate, this log file lists treatment/control pileup values, *p*- and *q*-values, and significance (`*`) for each interval. Here is part of the log file corresponding to the beginning of `peak_254` (Fig. 1):
```
chr1    1565338    1565605    0.000000    0.190111    0.000000    0.000000
chr1    1565605    1565608    1.000000    0.190111    2.284427    0.621277
chr1    1565608    1565618    2.000000    0.190111    4.568854    2.318294    *
chr1    1565618    1565647    3.000000    0.190111    6.853281    4.313010    *
chr1    1565647    1565654    2.000000    0.190111    4.568854    2.318294    *
chr1    1565654    1565720    1.000000    0.190111    2.284427    0.621277
chr1    1565720    1565733    2.000000    0.190111    4.568854    2.318294    *
```
* With multiple replicates, this log file lists *p*-values of each replicate, combined *p*-value, *q*-value, and significance for each interval.
* Note that this file (as well as the `-k` file, below) is called "bedgraph-ish" because it contains multiple `dataValue` fields, which isn't strictly allowed in the [bedGraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html).  However, a simple application of `awk` can produce the desired bedgraph files for visualization purposes (see this [awk reference](http://kirste.userpage.fu-berlin.de/chemnet/use/info/gawk/gawk_7.html#SEC57) for a guide to printing specific fields of input records).
<br><br>

```
  -k  <file>       Output bedgraph-ish file for pileups and p-values
```
* For each replicate, sequentially, this file lists a header line (`# treatment file: <name>; control file: <name>`), followed by treatment/control pileups and a *p*-value for each interval. This is the way to examine pileup values with multiple replicates, since the `-f` file will not supply them.
<br><br>

```
  -b  <file>       Output BED file for reads/fragments/intervals
```
* This is an unsorted [BED file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) of the reads/fragments/intervals analyzed. The 4th column gives the read name, number of alignments, 'T'reatment or 'C'ontrol, and sample number (0-based), e.g. `SRR5427885.57_2_T_0`.
<br><br>

## Filtering options <a name="filter"></a>

```
  -e  <arg>        Comma-separated list of chromosomes to ignore
```
* All alignments to the given list of chromosomes (reference sequences) will be ignored.  The alignments' lengths will not factor into the [total sequence information calculation](#pileup), nor to the average fragment length calculation (`-x`), and the alignments will not be printed to the `-b` file.
* For reads/fragments with multiple alignments, the scores of alignments to `-e` chromosomes *will* be considered for comparison purposes.
* The lengths of the `-e` chromosomes will be subtracted from the [total genome length](#genomelen) calculated by the program.
<br><br>

```
  -E  <file>       Input BED file(s) of genomic regions to ignore
```
* All alignments, or portions of alignments, that lie within the given genomic regions will be ignored.  The alignments' lengths (within an ignored region) will not factor into the [total sequence information calculation](#pileup).  However, the full fragment lengths *will* be counted for the average fragment length calculation (`-x`), and the full fragments *will* be listed in the `-b` file.
* The regions will affect peak calls, such that no peak may extend into or around an excluded region.
* In the output log files (`-f`, `-k`), excluded regions will have treatment/control pileup values of `0.0` and *p*-/*q*-values of `NA`.
* Multiple BED files can be specified, comma-separated (or space-separated, in quotes).  Overlapping BED intervals will be merged appropriately.
* The regions' lengths will be subtracted from the [total genome length](#genomelen) calculated by the program.
* Genomic regions to which reads typically do not align uniquely can be specified, but one should consider taking advantage of Genrich's ability to [analyze multimapping reads](#multimap).
* The accessory script [`findNs.py`](https://github.com/jsh58/Genrich/blob/master/findNs.py) will produce a BED file of 'N' homopolymers in a fasta file (e.g. a reference genome).  Its output should be given to Genrich via `-E`.
<br><br>

```
  -m  <int>        Minimum MAPQ to keep an alignment (def. 0)
```
* All alignments with `MAPQ` less than the given value will be eliminated.  This is equivalent to filtering with `samtools view -q <int>`.
* This option should not be used if the SAM/BAM lists [multiple alignments](#multimap) for some reads/fragments.  Instead, filtering should be accomplished via `-s <float>`, below.
<br><br>

<a name="sparam"></a>
```
  -s  <float>      Keep sec alns with AS >= bestAS - <float> (def. 0)
```
* Genrich analyzes all secondary alignments, but, by default, it keeps only the alignments whose [scores](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar) are equal to the best score for the read/fragment.  Setting a value such as `-s 20` will cause Genrich also to keep secondary alignments whose scores are within 20 of the best.
* The SAM/BAM should have alignment scores under the extra field `AS`.  If not, all alignments will be considered equivalent.
* Each of the `N` alignments for a read/fragment is counted as `1/N` for the pileup.
* As described [above](#multimap), a maximum of 10 alignments per read is analyzed.  Reads with more than 10 alignments within the `-s` threshold will be subsampled based on the best alignment scores; in the case of ties, alignments appearing first in the SAM/BAM are favored.
* The alignment score for a fragment (pair of reads) is equal to the sum of the reads' individual scores.
* Properly paired alignments take precedence over singleton alignments, regardless of the alignment scores.
<br><br>

### Unpaired alignments <a name="unpaired"></a>

By default, Genrich analyzes only properly paired alignments and infers the full fragments as spanning between the 5' ends of the two alignments (Fig. 2).  It does not analyze unpaired ("singleton") alignments unless one of three options is selected:
```
  -y               Keep unpaired alignments (def. false)
  -w  <int>        Keep unpaired alns, lengths changed to <int>
  -x               Keep unpaired alns, lengths changed to paired avg
```
* `-y`: unpaired alignments will be kept, just as they appear in the SAM/BAM
* `-w <int>`: unpaired alignments will be kept, with their lengths changed to the given value (from their 5' ends)
* `-x`: unpaired alignments will be kept, with their lengths changed to the average length of properly paired alignments (excluding those aligning to skipped chromosomes [`-e`])

<figure>
  <img src="figures/figure2.png" alt="Alignment analysis" width="700">
  <figcaption><strong>Figure 2.  Analysis of alignments by Genrich.</strong>  The alignment file <code>example.bam</code> has both properly paired alignments (top left) and unpaired "singleton" alignments (top right).  By default, Genrich infers the full fragments from the paired alignments and discards the unpaired alignments.  Unpaired alignments can be kept via <code>-y</code>, <code>-w &lt;int&gt;</code>, or <code>-x</code>, as described above.</figcaption>
</figure>
<br><br>

## ATAC-seq mode <a name="atacseq"></a>

[ATAC-seq](https://informatics.fas.harvard.edu/atac-seq-guidelines.html#overview) is a method for assessing genomic regions of open chromatin.  Since only the ends of the DNA fragments indicate where the transposase enzyme was able to insert into the chromatin, it may not be optimal to interpret alignments as shown above (Fig. 2).  Genrich has an alternative analysis mode for ATAC-seq in which it will create intervals centered on transposase cut sites (Fig. 3).

```
  -j               Use ATAC-seq mode (def. false)
  -d  <int>        Expand cut sites to <int> bp (def. 100)
```

<figure>
  <img src="figures/figure3.png" alt="ATAC-seq mode" width="700">
  <figcaption><strong>Figure 3.  ATAC-seq mode of Genrich.</strong>  Genrich analyzes intervals centered on cut sites (both ends of full fragments, as well as the 5' ends of unpaired alignments if <code>-y</code> is set).  The lengths of the intervals can be changed from the default of <code>-d 100</code>.</figcaption>
</figure>
<br><br>

Unpaired alignments can be analyzed with `-y`, though only one interval, centered on the read's 5' end, will be inferred.  Both `-w <int>` and `-x` are equivalent to `-y` in ATAC-seq mode.

The remainder of the peak-calling process (calculating pileups and significance values) is identical to the [default analysis mode](#method).  Note that the interval lengths (*not* the fragment lengths) are used to sum the total sequence information for the calculation of [control/background pileup values](#pileup).
<br><br>


## Peak-calling parameters<a name="peakcalling"></a>

```
  -q  <float>      Maximum q-value (FDR-adjusted p-value; def. 0.05)
  -p  <float>      Maximum p-value (overrides -q if set)
```
* This is the threshold below which a base is considered significantly enriched in the treatment sample(s) vs. the control/background.  Note that the significance value will be automatically converted to a -log<sub>10</sub> scale.
* When `-p` is selected, *q*-values will not be calculated (reported as -1).
<br><br>

```
  -a  <float>      Minimum AUC for a peak (def. 20.0)
```
* The default peak-calling method requires that, for a peak to be called, the total significance of the region must exceed a minimum value. The total significance is calculated as the sum of the -log(*q*) values above the `-q` threshold over the length of the region (i.e. the area under the -log(*q*) "curve").
* If a `-p` threshold is specified, the area under the -log(*p*) curve is calculated.
<br><br>

```
  -l  <int>        Minimum length of a peak (overrides -a if set)
```
* This option overrides the default peak-calling method (`-a`) and instead requires that peaks have at least the specified minimum length.  Any potential peaks whose lengths are below that threshold are eliminated, regardless of the significance.
<br><br>

```
  -g  <int>        Maximum distance between signif. sites (def. 100)
```
* This parameter sets the maximum distance between sites that achieve significance in order for them to be linked together into the same potential peak.  It applies in either peak-calling mode (`-a` and `-l`).
<br><br>


## Miscellaneous <a name="misc"></a>

```
  -z               Option to gzip-compress output(s)
```
* When selected, all output files will be gzip-compressed.
<br><br>

```
  -n  <int>        Number of threads to use (def. 1)
```
(not currently implemented)
<br><br>

Other options:
```
  -v/--verbose     Option to print status updates/counts to stderr
  -h/--help        Print the usage message and exit
  -V/--version     Print the version and exit
```
* Here is the verbose output for the sample processed above (Fig. 1):
```
$ ./Genrich  -t SRR5427884.bam  -o SRR5427884.narrowPeak  -f SRR5427884.log  -v
Processing treatment file #1: SRR5427884.bam
  BAM records analyzed:    6633684
    Paired alignments:     6633684
    Unpaired alignments:         0
  Fragments analyzed:      3316842
    Full fragments:        3316842
      (avg. length: 174.0bp)
- control file #1 not provided -
Peak-calling parameters:
  Genome length: 3036303846bp
  Significance threshold: -log(q) > 1.301
  Min. AUC: 20.000
  Max. gap between sites: 100bp
Peaks identified: 416514
```

### Warning messages <a name="warning"></a>

In verbose mode, Genrich may print one or more warnings to `stderr`:
* `Read N prevented from extending below 0 on <chrom>`: This may occur due to extending unpaired alignments (`-w <int>`, `-x`) or in ATAC-seq mode (`-j`).
* `Read N prevented from extending past <int> on <chrom>`: This also may occur due to extending unpaired alignments (`-w <int>`, `-x`) or in ATAC-seq mode (`-j`).
* `Large scaling may mask true signal`: This is printed if the [scaling factor](#pileup) for the control pileup is greater than 5.
* `BED interval ignored - located off end of reference`: An excluded BED interval (`-E`) whose start coordinate is past the end of the reference sequence will be ignored.  One should ensure that the genome version that produced the BED intervals matches that of the SAM/BAM.
* `BED interval extends past end of ref. - edited to <loc>`: An excluded BED interval (`-E`) whose end coordinate is past the end of the reference sequence will be adjusted as indicated.  Again, one should ensure that the genome version that produced the BED intervals matches that of the SAM/BAM.
* `No paired alignments to calculate avg frag length -- Printing singletons "as is"`: When there are *no* properly paired alignments and the [`-x` average extension option](#unpaired) is selected, the unpaired alignments will be printed as they appear in the SAM/BAM.
* `Read N, alignment at <loc> skipped due to overflow`: The maximum difference in pileup values from one genomic position to the next is +32767, and additional reads will be skipped due to this limitation.
* `Read N, alignment at <loc> skipped due to underflow`: The minimum difference in pileup values from one genomic position to the next is -32768, and additional reads will be skipped due to this limitation.
* `Read N has more than 128 alignments`: As described [above](#multimap), a maximum of 10 alignments per read is analyzed.  In addition, if a read has more than 128 alignments in the SAM/BAM, only the first 128 are considered.  Removing PCR duplicates may help reduce this issue.


## Contact <a name="contact"></a>

Genrich

Copyright &copy; 2018  John M. Gaspar (jsh58@wildcats.unh.edu)

