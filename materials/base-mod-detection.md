---
layout: page
element: notes
title: Base modification detection
---

In this session, we will perform reference-anchored modification calling
on the dataset we have been using thus far — yeast DNA where some thymidines have been
substituted by BrdU.
We will use the software program `DNAscent` on the files already generated
by us from previous sessions (fast5, bam, sequencing summary file).
We will then convert the output of DNAscent into the standard mod bam format for storing
modification information.
Along the way, we will also learn the mod bam file format.
We will be performing the steps highlighted with an asterisk in the pipeline figure below.

![Reference-unanchored pipeline with modification-calling highlighted](ref_unanc_workflow_modcall.png)

## Running DNAscent to produce modification calls along each sequenced DNA strand

### Preparations to call modifications

We first need to download a plugin to help DNAscent read the later versions of fast5 files.

```bash
mkdir -p ~/nanomod_course_references # any suitable folder will do here.
cd ~/nanomod_course_references
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -xvzf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
export HDF5_PLUGIN_PATH=$(pwd)/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
ls $HDF5_PLUGIN_PATH
# libvbz_hdf_plugin.so
rm ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz # if you see the output above, then cleanup by removing the tarball.
```

We now need to create an index needed by DNAscent — a plain text with two columns: read id
and the fast5 file with the corresponding time course of nanopore current.

```bash
DNAscent index -f ~/nanomod_course_data/carolin_nmeth_18/60 \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/index.dnascent \
  -s ~/nanomod_course_outputs/carolin_nmeth_18/sequencing_summary.txt
```

### Call modifications

We run the reference-anchored modification caller `DNAscent detect` here.
The program takes alignment information from the bam file, a reference genome,
and the raw nanopore currents as inputs and outputs the probability of modification
per thymidine per sequenced strand in a `.detect` format,
which is very similar to the TSV (tab-separated values) format.
The program uses a machine-learning approach where model parameters learned previously
represent differences in nanopore current characteristics between thymidine and BrdU.

The additional options below direct the program to use
16 computational threads and to reject alignments below a mapping quality of 20 and
a genome-mapped length of 1000 bases. 
By omitting the `--GPU` input parameter, we are running DNAscent in a slow, CPU-only mode
as the virtual machines used in the course do not have GPUs to lower costs.

```bash
DNAscent detect -b ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.onlyPrim.bam \
  -r ~/nanomod_course_references/sacCer3.fa \
 -i ~/nanomod_course_outputs/carolin_nmeth_18/index.dnascent \
 -o ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect \
 -t 16 \
 -q 20 -l 1000
```

## Call replication dynamics with DNAscent forkSense using single-molecule modification densities

<!-- TODO: Explain that we learn about features like forks in a previous session -->
<!-- TODO: Explain that forkSense is not strictly speaking needed for a general DNascent experiment -->
Run DNAscent forkSense

```bash
DNAscent forkSense -d ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect\
 -o ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.forkSense \
 -t 16 \
 --markOrigins --markTerminations --markForks
```

## Conversion of DNAscent detect into the modBAM format

We convert modification data from the specific `.detect` format used by DNAscent to
the more widely used and standardized `.mod.bam` or mod BAM format using the package
DNAscentTools developed by the Nieduszynski group.
Mod BAM files are just BAM files with two additional tags per line that store
modification information. We will discuss their structure later on in this session.

```bash
cd ~/nanomod_course_scripts/DNAscentTools
< ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect \
    python convert_detect_to_modBAM.py \
      --op ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.bam --tag T

# sort and index the BAM files
samtools sort -o ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.sorted.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.bam
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.sorted.bam
```

## Inspect modification data by conversion from modBAM to TSV format using modkit

At the end of our modification pipeline, all we want is a simple table with three columns
of read_id, coordinate, and modification_density.
An example with some artificial data is shown below.

```text
# marking thymidine modifications
read_id coordinate modification_density
f48c6a85-db3c-445f-865b-4bb876bd4a18 1000 0.1
f48c6a85-db3c-445f-865b-4bb876bd4a18 1004 0.4
f48c6a85-db3c-445f-865b-4bb876bd4a18 1014 0.9
...
5ddb8919-89de-4ffd-b052-423e53cff109 10210 0.2
5ddb8919-89de-4ffd-b052-423e53cff109 10224 0.3
5ddb8919-89de-4ffd-b052-423e53cff109 10249 0.88
```

We can use the `modkit extract` program to convert data in our mod BAM file into a tabular format.
Please run the code below and inspect the output `tsv` file.
You should see columns such as `read_id`, `forward_read_position`,
`ref_position`, `chrom`, `mod_qual` etc.

```bash
cd ~/nanomod_course_outputs/carolin_nmeth_18/
modkit extract dnascent.detect.mod.sorted.bam dnascent.detect.mod.sorted.bam.tsv
```

## Discussion of the modBAM file format

We conclude this session with a [discussion]({{ site.baseurl }}/materials/mod-bam-format)
 of the mod BAM file format.