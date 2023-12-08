---
layout: page
element: notes
title: Sequence alignment and pycoQC
---

TBD sequence alignment and pycoQC.

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