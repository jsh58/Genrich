# Genrich: detecting sites of genomic enrichment

## Table of Contents
* [Introduction](#intro)
  * [Quick start](#quick)
  * [Software compilation](#compile)
  * [Usage message](#usage)
* [Peak-calling method](#method)
* [I/O files and options](#files)
  * [Required files](#required)
  * [Optional files](#optional)
* [Miscellaneous](#misc)
* [Contact](#contact)
<br><br>

## Introduction <a name="intro"></a>

Genrich is a peak-caller for genomic enrichment assays (e.g. ChIP-seq, ATAC-seq).  It analyzes alignment files generated following the assay and produces a peak file.
<br><br>

### Quick start <a name="quick"></a>

Given:
* `sample.bam` (alignment file)
* `Genrich` (compiled as described [below](#compile))
<br><br>

To produce a file listing regions of genomic enrichment:
```
$ ./Genrich  -t sample.bam  -o sample.narrowPeak
```
<br>

### Software compilation <a name="compile"></a>

The software can be downloaded from [GitHub](https://github.com/jsh58/Genrich).

A Makefile is provided for compilation with [GCC](https://gcc.gnu.org/releases.html), and [zlib](http://zlib.net) is also required.  The program has been tested after compilation with GCC 5.4.0 and zlib 1.2.8.

To compile, run `make` in the folder in which the software was downloaded.  The executable `Genrich` should be produced.
<br><br>


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
<br>

## Peak-calling method <a name="method"></a>

(coming soon)


## I/O files and options <a name="files"></a>

### Required files <a name="required"></a>

```
  -t  <file>       Input SAM/BAM file(s) for treatment sample(s)
```
Genrich analyzes alignment files in [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).  SAM files must have a header.

SAM/BAM files for multiple replicates can be specified, comma-separated (or space-separated, in quotes).  Multiple SAM/BAM files for a single replicate should be combined in advance via `samtools merge`.

The SAM/BAM files should be name sorted (via `samtools sort -n`).  As of [Version 0.3](https://github.com/jsh58/Genrich/releases/tag/v0.3), unsorted SAM/BAM files are allowed, but this is likely to change.
<br><br>

```
  -o  <file>       Output peak file (in ENCODE narrowPeak format)
```
As indicated, the output file is in [ENCODE narrowPeak format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12).  Here are additional details of the fields:
<table>
  <tr>
    <td align="center">4. name</td>
    <td><code>peak_N</code>, where <code>N</code> is the 0-based count</td>
  </tr>
  <tr>
    <td align="center">5. score</td>
    <td>10 * qValue (or 10 * pValue if <code>-p</code> is set), rounded to the nearest int (max. 1000)</td>
  </tr>
  <tr>
    <td align="center">7. signalValue</td>
    <td>Total area under the curve, in default peak-calling mode.  In minimum-length peak-calling mode (<code>-l</code>), the summit fold-enrichment (treatment / control) with one replicate, or the pValue/qValue with multiple replicates.</td>
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
    <td>Summit position: the midpoint of the interval reaching the highest significance (the longest interval in case of ties)</td>
  </tr>
</table>
<br><br>

### Optional files <a name="optional"></a>

```
  -c  <file>       Input SAM/BAM file(s) for control sample(s)
                     (matched with -t files; 'null' if missing)
```
Alignment files for control samples can be specified.  As indicated, they should be matched with treatment files.
<br><br>

```
  -f  <file>       Output bedgraph-ish file for p/q values
```
When Genrich analyzes a single replicate, this log file lists treatment/control pileup values, *p*- and *q*-values, and significance (`*`) for each interval.  With multiple replicates, it lists *p*-values of each replicate, combined *p*-value, *q*-value, and significance.

Note that this file (as well as the `-k` file, below) is called "bedgraph-ish" because it contains multiple `dataValue` fields, which isn't strictly allowed in the [bedGraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html).  However, a simple application of `awk` can produce the desired bedgraph files for visualization purposes.
<br><br>

```
  -k  <file>       Output bedgraph-ish file for pileups and p-values
```
For each replicate, sequentially, this file lists a header line (`# treatment file: <name>; control file: <name>`), followed by treatment/control pileups and a *p*-value for each interval. This is the way to examine pileup values with multiple replicates, since the `-f` file will not supply them.
<br><br>

```
  -b  <file>       Output BED file for reads/fragments/intervals
```
This is an unsorted BED file of the reads/fragments/intervals analyzed. The 4th column gives the read name, number of alignments, 'T'reatment or 'C'ontrol, and sample number (0-based), e.g. `SRR5427885.57_2_T_0`.
<br><br>

## Miscellaneous <a name="misc"></a>

```
  -z               Option to gzip-compress output(s)
```
When selected, all output files will be gzip-compressed.
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
<br><br>

## Contact <a name="contact"></a>

Genrich

Copyright &copy; 2018  John M. Gaspar (jsh58@wildcats.unh.edu)

