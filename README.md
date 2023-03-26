# **Introduction**

This is the custom code used in the following paper: **Epigenetic therapy activates TE-chimeric transcripts to provide additional source of antigens in glioblastoma**

Once the paper is published, the citation information will be here: X

Of note, the transcript analysis pipeline used in this paper can be found here: https://github.com/twlab/TEProf2Paper

Furthermore, the script used to perform cap filtering for nanocage data can be found here: https://github.com/nakul2234/cageCapfilter

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
