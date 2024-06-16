# **Introduction**

This is the custom code used in the following paper: **Epigenetic therapy potentiates transposable element transcription to create tumor-enriched antigens in glioblastoma cells**

Once the paper is published, the citation information will be here: X

Of note, the transcript analysis pipeline used in this paper can be found here: https://github.com/twlab/TEProf2Paper

Furthermore, the script used to perform cap filtering for nanocage data can be found here: https://github.com/nakul2234/cageCapfilter

The script for peak calling using nanoCAGE data and transcript profiling using long-read data can be found here: https://github.com/twlab/LRCAGE

All scripts defined below can be found in the `bin` folder.

# **Scripts and Descriptions**

## TITEA_Analysis.Rmd

R Notebook with analysis of mRNA, nanoCAGE, long read sequencing, and mass spectrometry data

## rmskAllFramesTranslate.py

Python script that translates repeatmasker fasta files in all frames to create a potential protein database. 

```
usage: rmskALLFramesTranslate.py [-h] [-m] [-e] [-l <LENGTH>] <FASTA FILE>

Translate FASTA file in all 6 frames. Optimized for use with Repeatmasker TE
files, but any FASTA can be tried. The FASTA fieldname orders will be
different.

positional arguments:
  <FASTA FILE>  The fasta file with DNA sequences of TEs that will be
                translated by this program. This script has been optimized to
                work with Repeatmasker Files.

optional arguments:
  -h, --help    show this help message and exit
  -m            If this flag is placed, then only those peptides starting with
                a Methionine will be considered
  -e            If this flag is placed, then peptides with no stop codon will
                be considered
  -l <LENGTH>   Minimum peptide length. Default: 1
```

## transcriptALLFramesTranslate.py

```
usage: transcriptALLFramesTranslate.py [-h] [-m] [-e] [-l <LENGTH>]
                                       <FASTA FILE>

Translate FASTA file in all 3 forward frames. Optimized for use with FASTA
files with a uniqid as the only field after >, but any FASTA can be tried.

positional arguments:
  <FASTA FILE>  The fasta file with DNA sequences that will be translated by
                this program.

optional arguments:
  -h, --help    show this help message and exit
  -m            If this flag is placed, then only those peptides starting with
                a Methionine will be considered
  -e            If this flag is placed, then peptides with no stop codon will
                be considered
  -l <LENGTH>   Minimum peptide length. Default: 1
```

## call_ATACpeak.sh
```
usage: call_ATACpeak.sh <ATAC BAM> <BLACKLIST BED>
Call ATAC peaks from ATAC bam file. To infer Tn5 insertion site for each read, +5 and -4 offset is applied to reads aligned to the forward and reverse strand, respectively.

<positional arguments>:
  <ATAC BAM>       The bam file of ATAC-seq data.
  <BLACKLIST BED>  The blacklist bed file that will filter out overlapping ATAC peaks.
```

## calc_ATACscalefactor.R
```
usage: calc_ATACscalefactor <ATAC SAMPLE SHEET CSV> <ATAC PEAK AT HOUSE KEEPING GENES BED> <SIZE FACTOR TXT>
Calculate scale factors for each ATAC data using consensus ATAK peaks at house keeping gene promoters. 

<positional arguments>:
  <ATAC SAMPLE SHEET CSV>                 Sample sheet CSV file of ATAC samples. This format is used for DiffBind package.
  <ATAC PEAK AT HOUSE KEEPING GENES BED>  The BED file of consensus ATAC peaks overlapping house keeping gene promoters.
  <SIZE FACTOR TXT>                       Scale factors calculated for each ATAC sample.
```


## pipe_align_ATAC.sh
```
usage: pipe_align_ATAC.sh <ATAC READ 1> <ATAC READ 2>
Trim adaptors in ATAC-seq reads, align to the Hg38 reference genome using bowtie2, deduplicate, extract usable reads, call peaks, and filter blacklist regions.

<positional arguments>:
  <ATAC READ 1> ATAC-seq read 1.
  <ATAC READ 2> ATAC-seq read 2.
```

## pipe_align_WGBS.sh
```
usage: pipe_align_WGBS.sh <WGBS READ 1> <WGBS READ 2> <OUTPUT DIR>
Trim adaptors in WGBS-seq reads, align to the Hg38 reference genome using bismark, deduplicate, and extract methylation calls.

<positional arguments>:
  <WGBS READ 1>     WGBS-seq read 1.
  <WGBS READ 2>     WGBS-seq read 2.
  <OUTPUT DIR> Output directory. 
```

## pipe_align_nanoCAGE.sh
```
usage: pipe_align_nanoCAGE.sh <nanoCAGE READ 1> <nanoCAGE READ 2>
Trim adaptors in nanoCAGE-seq reads, align to the Hg38 reference genome using STAR.

<positional arguments>:
  <nanoCAGE READ 1>     nanoCAGE-seq read 1.
  <nanoCAGE READ 2>     nanoCAGE-seq read 2.
```

## pipe_align_RNA.sh
```
usage: pipe_align_RNA.sh <RNA READ 1> <RNA READ 2>
Trim adaptors in RNA-seq reads, align to the Hg38 reference genome using STAR.

<positional arguments>:
  <RNA READ 1>     RNA-seq read 1.
  <RNA READ 2>     RNA-seq read 2.
```

## pipe_assemble_RNA.sh
```
usage: pipe_assemble_RNA.sh <RNA BAM>
Assemble transcripts with stringtie using RNA-seq BAM file as an input (only uniquely mapped reads are used).

<positional arguments>:
  <RNA BAM> RNA BAM file.
```




