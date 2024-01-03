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

## Running DNAscent to produce modification calls along each sequenced DNA strand

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

Run DNAscent index

```bash
DNAscent index -f ~/nanomod_course_data/carolin_nmeth_18/60 \
  -o ~/nanomod_course_outputs/carolin_nmeth_18/index.dnascent \
  -s ~/nanomod_course_outputs/carolin_nmeth_18/sequencing_summary.txt
```

Run DNAscent detect

```bash
DNAscent detect -b ~/nanomod_course_outputs/carolin_nmeth_18/aligned_reads.sorted.onlyPrim.bam \
  -r ~/nanomod_course_references/sacCer3.fa \
 -i ~/nanomod_course_outputs/carolin_nmeth_18/index.dnascent \
 -o ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect \
 -t 16 \
 -q 20 -l 1000
```

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

Make mod bam file

```bash
cd ~/nanomod_course_scripts/DNAscentTools
< ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect \
    python convert_detect_to_modBAM.py \
      --op ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.bam --tag T

samtools sort -o ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.sorted.bam \
  ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.bam
samtools index ~/nanomod_course_outputs/carolin_nmeth_18/dnascent.detect.mod.sorted.bam
```

## Discussion of the modBAM file format

<!-- TODO -->

## Conversion from modBAM to TSV format using modkit

<!-- TODO -->

Get tabular output from mod bam file

```bash
cd ~/nanomod_course_outputs/carolin_nmeth_18/
modkit extract dnascent.detect.mod.sorted.bam dnascent.detect.mod.sorted.bam.tsv
```