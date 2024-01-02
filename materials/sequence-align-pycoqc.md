---
layout: page
element: notes
title: Sequence alignment and pycoQC
---

In this session, we will convert raw currents per DNA strand recorded by nanopore devices
into DNA sequences (basecalling), find the best-fit location for each strand on a
reference genome (alignment), and assess their quality (quality control).
These are the first steps in any nanopore experiment and are followed by specialized
steps such as genome assembly or modification calling.

## Basecalling: converting nanopore currents into DNA sequences

Nanopore devices ship with basecalling software that converts nanopore currents into DNA
sequences. Basecallers use pre-recorded model parameters that contain information about the
characteristic currents produced by different DNA k-mers as they translocate through the nanopore.
These programs segment the current corresponding to one DNA strand, assign a k-mer per segment,
and stitch k-mers together into one sequence. We will use two basecallers: `guppy` in this session
and the more-recent `dorado` in a later session.

Basecallers can be run on the nanopore device during sequencing and/or can be run later on a computer.
For most purposes, on-device basecalling is sufficient. If a higher accuracy is desired or
better model parameters become available, then data can be re-basecalled on a computer.
Accuracy is specified through input parameters and higher accuracy leads to longer runtimes.

Basecallers generally output DNA sequences with the four canonical DNA bases;
modification calling is a separate step and is bolted on to the output of basecallers.
We will deal with modification calling in a later session.

### Input and output file formats used by basecallers

Basecallers use nanopore currents and model parameters as inputs and produce one DNA sequence
per strand as output. Here we will look at file formats used to store nanopore currents
and DNA sequences. We will not discuss how model parameters are stored.

#### Fast5/pod5 files contain nanopore currents

Historically, raw nanopore currents were stored in `.fast5` files.
Each read (DNA strand) is assigned a unique 36 character string (read id)
e.g. `acde070d-8c4c-4f0d-9d8a-142843c10333` and a current time course as
well as some metadata such as the label of the pore that sequenced this strand.
Fast5 files can contain one read per file or multiple reads.
Recently, fast5 files were replaced by the more efficient `.pod5` files.
We will use fast5 files in this session and pod5 files in a later session.
These are not plain text files and cannot be inspected directly on the command line.

#### Fasta/fastq/bam files contain DNA sequences

Fasta files are in plain text and store many sequences per file.
They are the simplest, with one line per read id followed by one
line of sequence. For readability, the sequence line may be
broken up into several lines, with all but the last line stipulated
to have the same width.

```text
>acde070d-8c4c-4f0d-9d8a-142843c10333
TTGTA
```

Fastq files are plain text files and contain more information,
with four lines per sequenced read, e.g.:

```text
@acde070d-8c4c-4f0d-9d8a-142843c10333 ch=139 start_time=2022-09-26T18:43:42Z
TTGTA
+
8,'&)
```

- The first line contains the read id and associated metadata in key-value pairs after the `@` symbol
- The second line is the sequence itself. Please note that the entire sequence is on one line unlike a fasta file.
- The third line has a `+` symbol and may contain other information.
- The fourth line contains quality information per called base represented as ASCII characters.
  For example, the fifth character `)` associated with the last `A` has an ASCII code of 41.
  Quality increases exponentially with the numeric code: the probability of a wrong basecall drops
  ten fold for an increment of 10 in the ASCII code.

BAM files are the most flexible, as they may contain one or all of the basecalling/alignment/modification
information per read. We will not discuss them here but later on in this session when we discuss alignment.

We will use the `.fastq` format in our pipeline below to store basecalling outputs,
and the `.fasta` format to store the reference genome for alignment purposes.

## Redistribute the material below.

Download the sacCer3 reference genome and make a fasta index.

```bash
mkdir -p ~/nanomod_course_references # any folder will do here.
cd ~/nanomod_course_references
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
samtools faidx sacCer3.fa
```

TBD: Downloading carolin's dataset and setting up just the `60/` fast5 folder 
at `~/nanomod_course_data/carolin_nmeth_18`

We do the basecalling using guppy

```bash
mkdir -p ~/nanomod_course_outputs/carolin_nmeth_18/

guppy_basecaller --input_path ~/nanomod_course_data/carolin_nmeth_18/60 \
    --save_path ~/nanomod_course_outputs/carolin_nmeth_18 \
    --config dna_r9.4.1_450bps_hac.cfg --recursive -q 10 --disable_pings \
    --num_callers 8 --cpu_threads_per_caller 2 --disable_qscore_filtering \
    --progress_stats_frequency 10
```

We stop the job after five minutes or so, after we've four fastq files.

Minimap2 alignment

```bash
minimap2 -L -a -x map-ont -t 8 \
  ~/nanomod_course_references/sacCer3.fa \
  ~/nanomod_course_outputs/carolin_nmeth_18/*.fastq.* >\
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sam
```

Samtools sort and index

```bash
samtools sort -@ 16 -T /tmp \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sam
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam
```

Run pycoQC

```bash
pycoQC \
    -f ~/nanomod_course_outputs/carolin_nmeth_18/sequencing_summary.txt \
    -a ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam \
    -o ~/nanomod_course_outputs/carolin_nmeth_18/pycoQC/analysis.html \
    -j ~/nanomod_course_outputs/carolin_nmeth_18/pycoQC/analysis.json \
    --quiet
```

Motivate that we retain only primary reads.

```bash
samtools view -Sb -F "$((256 + 2048))" \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.onlyPrim.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam;
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam;
```