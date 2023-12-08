---
layout: page
element: notes
title: Comparison of methods for single molecule base modification detection
---

Section in progress.

## Download input data

First, we download the fast5 files used as input in this section and convert them to pod5 format.

```bash
mkdir -p ~/nanomod_course_data   # this can be any folder where input data is stored
cd ~/nanomod_course_data
mkdir -p 20220503_1708_1H_PAM63103_c3de1ce7
cd 20220503_1708_1H_PAM63103_c3de1ce7
aws s3 cp --no-sign-request s3://ont-open-data/cliveome_kit14_2022.05/cfdna/flowcells/20220503_1708_1H_PAM63103_c3de1ce7/fast5/PAM63103_0cd67e80_10000.fast5 .

mkdir -p ~/nanomod_course_outputs  # this can be any folder where outputs are stored
cd ~/nanomod_course_outputs
mkdir -p 20220503_1708_1H_PAM63103_c3de1ce7
cd 20220503_1708_1H_PAM63103_c3de1ce7
pod5 convert fast5 ~/nanomod_course_data/20220503_1708_1H_PAM63103_c3de1ce7/PAM63103_0cd67e80_10000.fast5 --output PAM63103_0cd67e80_10000.pod5
```

We can inspect the pod5 file to see how many reads are on it.

```bash
pod5 inspect summary PAM63103_0cd67e80_10000.pod5
```

There are 4000 reads here. The output should look like the following.

```text
File version in memory 0.3.2, read table version 3.
File version on disk 0.3.2.
File uses VBZ compression.
Batch 1, 1000 reads
Batch 2, 1000 reads
Batch 3, 1000 reads
Batch 4, 1000 reads
Found 4 batches, 4000 reads
```

## Call modifications

We need to download the models `dorado` uses for basecalling.
First, we make a suitable directory for these models

```bash
mkdir -p ~/nanomod_course_references  # this can be any folder where reference information is stored
cd ~/nanomod_course_references        
mkdir -p dorado_models
```

Then, we download the models.

```bash
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v3.5.2 --directory ./dorado_models
dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.1.0 --directory ./dorado_models
```

We make a directory to store our basecalls.

```bash
mkdir -p ~/nanomod_course_outputs  # this can be any folder where outputs are stored
cd ~/nanomod_course_outputs
mkdir -p 20220503_1708_1H_PAM63103_c3de1ce7
cd 20220503_1708_1H_PAM63103_c3de1ce7
```

Now, we basecall 10 reads (`-n 10`) from the cliveome fast5 file using `dorado` without modification calling. NOTE: the `-b 10 -c 1000` are internal `dorado` parameters which we've chosen to fit our computer.

```bash
dorado basecaller \
    ~/nanomod_course_references/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v3.5.2 \
    ./ \
    --verbose -n 10 -x cpu -b 10 -c 1000 | \
    samtools view --threads 8 -O BAM -o ./PAM63103_0cd67e80_10000.only_10_reads.bam
```

We can basecall 10 reads with 5mC methylation.

```bash
dorado basecaller \
    ~/nanomod_course_references/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v3.5.2 \
    ./ \
    --verbose -n 10 -x cpu -b 10 -c 1000 --modified-bases 5mCG | \
    samtools view --threads 8 -O BAM -o ./PAM63103_0cd67e80_10000.only_10_reads.5mCG.bam
```

We can basecall 10 reads with both 5hmC and 5mC methylation.

```bash
dorado basecaller \
    ~/nanomod_course_references/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v4.1.0 \
     ./ \
    --verbose -n 10 -x cpu -b 10 -c 1000 --modified-bases 5mCG_5hmCG | \
    samtools view --threads 8 -O BAM -o ./PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam
```

We can view the output with `modkit`.

```bash
samtools index PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam
modkit extract --force PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam.tsv
# --force makes modkit overwrite output files files if they exist
wc -l PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam.tsv
# 85
head -n 10 PAM63103_0cd67e80_10000.only_10_reads.5mCG_5hmCG.bam.tsv
```

This will produce the following output

```text 
read_id forward_read_position   ref_position    chrom   mod_strand      ref_strand      ref_mod_strand  fw_soft_clipped_start   fw_soft_clipped_end     read_length     mod_qual        mod_code        base_qual r
ef_kmer query_kmer      canonical_base  modified_primary_base   inferred
0058c812-6b1e-459d-a82e-c01bba1e8ff5    30      -1      .       +       .       .       0       0       396     0.7089844       h       5       .       TACGT   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    30      -1      .       +       .       .       0       0       396     0.001953125     m       5       .       TACGT   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    142     -1      .       +       .       .       0       0       396     0.005859375     h       42      .       GGCGT   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    142     -1      .       +       .       .       0       0       396     0.9902344       m       42      .       GGCGT   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    158     -1      .       +       .       .       0       0       396     0.087890625     h       47      .       CCCGG   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    158     -1      .       +       .       .       0       0       396     0.9082031       m       47      .       CCCGG   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    292     -1      .       +       .       .       0       0       396     0.013671875     h       12      .       TTCGC   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    292     -1      .       +       .       .       0       0       396     0.009765625     m       12      .       TTCGC   C       C       false
0058c812-6b1e-459d-a82e-c01bba1e8ff5    300     -1      .       +       .       .       0       0       396     0.029296875     h       21      .       ATCGG   C       C       false
```