# **Introduction**

This is the custom code used in the following paper: **Epigenetic therapy activates TE-chimeric transcripts to provide additional source of antigens in glioblastoma**

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
  <SIZE FACTOR TXT>                       Scale factors calculated for each ATAC sample
```

