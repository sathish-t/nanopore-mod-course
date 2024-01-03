---
layout: page
element: notes
title: Sequence alignment and pycoQC
---

In this session, we will convert raw currents per DNA strand recorded by nanopore devices
into DNA sequences (basecalling), find the best-fit location for each strand on a
reference genome (alignment/genome mapping), and assess their quality (quality control).
Basecalling and quality control are the first steps in any nanopore experiment,
and are followed by more specialized steps like genome assembly or modification calling.

<!-- TODO: pipeline figure -->

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

In this section, we will first look at the input and output file formats used by basecallers,
and then run the basecalling commands on nanopore data from yeast.

### Input and output file formats used by basecallers

Basecallers use nanopore currents and model parameters as inputs and produce one DNA sequence
per strand as output. Here we will look at different file formats used to store nanopore currents
and DNA sequences. You can choose the file formats you use in your experiments.
We will not discuss how model parameters are stored.

#### Fast5/pod5 files contain nanopore currents

Historically, raw nanopore currents were stored in `.fast5` files.
In them, each read (DNA strand) has a unique 36 character string (read id)
e.g. `acde070d-8c4c-4f0d-9d8a-142843c10333` and a current time course as
well as some metadata such as the label of the pore that sequenced this strand.
Fast5 files can contain one read per file or multiple reads.
Recently, fast5 files were replaced by the more efficient `.pod5` files.
We will use fast5 files in this session and pod5 files in a later session.
These are not plain text files and cannot be inspected directly on the command line.

#### Fasta/fastq/bam files contain DNA sequences

One can use several file formats to store DNA sequence data.
Fasta files can only contain read ids and sequences, whereas fastq files
can include additional information like the quality of basecalling per base.
BAM files are the most flexible, as they may contain one or all of the basecalling/alignment/modification
information per read. We use the `.fastq` format in our pipeline below to store basecalling outputs,
and the `.fasta` format to store the reference genome for alignment purposes.
We will not discuss BAM files here but later on in this session when we discuss alignment.

Fasta files are in plain text and store many sequences per file.
They are the simplest, with one line per read id followed by one
line of sequence. For readability, the sequence line may be
broken up into several lines, with all but the last line stipulated
to have the same width. Fasta files often have an associated index
ending in `.fai` (e.g.: `something.fa.fai` accompanying `something.fa`)
which enables fast retrieval of sequence data. An example line from
a fasta file is shown below.

```text
>acde070d-8c4c-4f0d-9d8a-142843c10333
ATCGA
```

Fasta files are commonly used to store reference genomes. For example,
let's say an organism had two chromosomes which are both very short.
The fasta file of its reference genome may look like

```text
>chrI
gatgcaaagcatcggcttttactcacgatccgacgacacaattcagcgacagggcactcc
caaatcacgcttgcgggaaatatcactctt
>chrII
tcctcacgaaattaactagaagacccagcatatccaacaagggatgagtgtcacttacgc
caatcctcgtctcgagcttatcggtttgta
```


Fastq files are plain text files and contain more information,
with four lines per sequenced read, e.g.:

```text
@acde070d-8c4c-4f0d-9d8a-142843c10333 ch=139
ATCGA
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

### Running the basecalling commands

We do the basecalling using guppy as shown below. The options specify input file paths,
output file paths, the model file of choice (.cfg), and details on how many reads
per fastq file and how many threads of execution. Normally, reads are split into pass
and fail bins depending on their quality, but we allow all reads through as our
reads contain modified bases which are expected to interfere with basecalling accuracy.

<!-- TODO: TBD: Downloading carolin's dataset and setting up just the `60/` fast5 folder 
at `~/nanomod_course_data/carolin_nmeth_18` -->

```bash
mkdir -p ~/nanomod_course_outputs/carolin_nmeth_18/

guppy_basecaller --input_path ~/nanomod_course_data/carolin_nmeth_18/60 \
    --save_path ~/nanomod_course_outputs/carolin_nmeth_18 \
    --config dna_r9.4.1_450bps_hac.cfg --recursive -q 10 --disable_pings \
    --num_callers 8 --cpu_threads_per_caller 2 --disable_qscore_filtering \
    --progress_stats_frequency 10
```

We stop the job after five minutes or so using `Ctrl+C` (`Command+.` for mac users),
after we've four fastq files.

## Alignment: locating sequenced DNA on a reference genome

Pre-existing reference genomes assembled by others are available from databases for
some organisms. Alignment/genome mapping is the process of matching each sequenced
DNA strand to locations along the reference. 
In a modification-calling experiment, alignment helps correlate DNA modification density with genomic
locations and/or genomic features, so that one can answer questions like 
'which part of the genome is highly modified?', 'are modification densities higher within genes?' etc.
We will use the aligner `minimap2`,
which takes experimental sequence data and a linear reference genome as input and
reports alignments in the BAM format as requested by us.
We will give an overview of alignment and the file formats used before running the commands.

In a modification-calling genomics pipeline,
alignment can be performed before or after modification calling. 
In our workflow today, we will use a program that implements reference-anchored modification
calling, so we need to do alignment before calling modifications.
The program, DNAscent, calls modifications on the reference sequence corresponding
to each read instead of the basecalled sequence.
Tomorrow, we will use a program that performs reference-independent modification calling
and performs alignment after if the user requests so.

### Aligners output multiple alignments per strand with alignment quality and per-base alignment operations

Alignment programs may output zero to several alignments per DNA sequence.
If multiple alignments are output, the ones other than the best alignment
are labelled as secondary or supplementary alignments.
We retain only primary alignments for modification-calling purposes.

It is unusual for a strand to match
perfectly with a location on the reference, so the aligner usually outputs
a quality score along with the region on the reference and
a series of operations that must be performed on the sequence to match it with the region.
For example, let's say an aligner has determined that out of all possible regions on the reference,
the read `acde070d-8c4c-4f0d-9d8a-142843c10333` with sequence `ATCGA`
(a read is also called a query by aligners) matches best with the region `chr1:10-15` of
the sequence `ATCGC`. The aligner will report match information like `4M1X` which means the
first four bases were a match but the last base was a mismatch.
We will not discuss details of how an aligner works or how to read these so-called CIGAR strings
that report alignment information in this course.

### Input and output file formats used by minimap2

Aligners can use multiple file formats for input and output.
In our pipeline, we will feed DNA sequences in the `fastq` format and use a linear
reference genome in the `fasta` format, both of which we've already discussed above.
We discuss the output BAM file format below.

### BAM/SAM file format is used to store alignment information

BAM files are binary files that store alignment information.
SAM files are BAM files but rendered in plain text, and have two sections:
- a header section containing overall information about the file/the workflow, 
where each line begins with an `@` symbol 
- an alignment section with one line per alignment, where each line
begins with the read id and columns are tab-separated.
Note that one read id can have multiple alignments.

The alignment line of the read `acde070d-8c4c-4f0d-9d8a-142843c10333` looks like 
the following in plain text (tabs have been replaced with spaces for readability).

```text
acde070d-8c4c-4f0d-9d8a-142843c10333 0 chr1 11 50 4M1X * 0 0 ATCGA 9</(; XR:i:975
```

We will not go through all the columns in detail. The important ones are:
- the first column is the read id
- the second column is a flag and contains several bits of information, including whether the alignment is secondary/supplementary.
- the third column is the contig on the reference genome
- the fourth column is the starting position of the mapping on the reference genome (in 1-based coordinates)
- the fifth column is the quality of mapping
- the sixth column is the CIGAR string
- the tenth column is the DNA sequence of the read

The first 11 columns are mandatory and are followed by optional tags with the TAG:TYPE:VALUE format.
Using the optional `ML` and `MM` tags, one can store modification data.
We will discuss BAM files with modification information, called modbam files, in a later session.
For more information about the columns and the optional tags, please consult the file format
specifications [here](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[here](https://samtools.github.io/hts-specs/SAMtags.pdf).

### (optional) BAM/SAM file format is versatile

BAM files are versatile and can store many types of data. The format is the output of choice
for many other types of software programs such as basecallers and modification callers as well.
Their outputs would deviate from the example line shown above by removing bits of data.
For e.g. one can forego alignment information altogether by setting the flag to 4 and
a few other columns to 0 or *, or forego basecalling information by setting the sequence to *.

### Samtools: a collection of programs that operate primarily on BAM files

Samtools is a collection of programs that perform a wide variety of tasks on
BAM files e.g. indexing them, subsetting them in a random manner or
by reads that map to a genomic region etc. Samtools can also convert
between different file formats e.g. BAM to SAM, FASTQ to BAM etc.
We will use several such commands in our workflow and
we will introduce them as and when we use them.
The commands generally have the syntax `samtools <command> <input_file>`.

### Running the alignment commands

We will first obtain the sacCer3 (_S. cerevisiae_) reference genome
and make a corresponding fasta index file.

```bash
mkdir -p ~/nanomod_course_references # any folder will do here.
cd ~/nanomod_course_references
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
samtools faidx sacCer3.fa
```

We perform alignment using `minimap2` with 8 threads, using the basecalled fastq
files and the reference genome we downloaded above. The `-x map-ont` parameter instructs
minimap2 to use alignment parameters optimized for nanopore data.

```bash
minimap2 -L -a -x map-ont -t 8 \
  ~/nanomod_course_references/sacCer3.fa \
  ~/nanomod_course_outputs/carolin_nmeth_18/*.fastq.* >\
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sam
```

Whenever BAM/SAM files are made, it is a good idea to sort and index them as many
tools require that BAM/SAM files are sorted and indexed.
We can use samtools to do this.

```bash
samtools sort -@ 16 -T /tmp \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sam
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam
```

### Inspecting alignments using genome browsers and samtools

<!-- TODO: complete this -->

## Quality control

Run pycoQC

```bash
pycoQC \
    -f ~/nanomod_course_outputs/carolin_nmeth_18/sequencing_summary.txt \
    -a ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam \
    -o ~/nanomod_course_outputs/carolin_nmeth_18/pycoQC/analysis.html \
    -j ~/nanomod_course_outputs/carolin_nmeth_18/pycoQC/analysis.json \
    --quiet
```

You can open the analysis.html file in your browser after the program has finished running.
You should see a webpage whose layout and figures, but not the actual details, are similar to
[this](https://a-slide.github.io/pycoQC/pycoQC/results/Guppy-2.3_basecall-1D_alignment-DNA.html).

<!-- TODO: talk about what you see when you run pycoQC -->

## Filter BAM file to include only primary reads

<!-- TODO: Motivate that we retain only primary reads, only short section needed here --> 

```bash
samtools view -Sb -F "$((256 + 2048))" \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.onlyPrim.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam;
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.bam;
```