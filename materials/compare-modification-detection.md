---
layout: page
element: notes
title: Alternative method for single molecule base modification detection
---

In this session, we will perform reference-independent modification calling
on a dataset of human DNA that naturally contains cytosines modified due to methylation.
We will use the basecaller `dorado` that can perform basecalling and modification
calling when we issue one command.
We will be using a subset of the 'Cliveome' [dataset](https://labs.epi2me.io/cliveome_5mc_cfdna_celldna/)
developed by Oxford NanoporeTech.
We have already downloaded the input data and copied it to the location `~/nanomod_course_data/human/`
in a [previous]({{ site.baseurl }}/data) step.

## Inspect input data

We can inspect the pod5 file to see how many reads are on it.

```bash
input_pod5=~/nanomod_course_data/human/PAM63974_pass_58881fec_60.twenty_random_reads.pod5
pod5 inspect summary $input_pod5
```

There are 20 reads here. The output should look like the following.

```text
File version in memory 0.3.2, read table version 3.
File version on disk 0.3.2.
File uses VBZ compression.
Batch 1, 20 reads
Found 1 batches, 20 reads
```

## Basecalling and modification calling

We need to download the models `dorado` uses for basecalling and modification calling.
We will put them in a suitable directory.

```bash
dorado_model_dir=~/nanomod_course_references/dorado_models
model_config=dna_r10.4.1_e8.2_400bps_hac@v3.5.2
mkdir -p $dorado_model_dir
dorado download --model $model_config --directory $dorado_model_dir
```

We make a directory to store our modification calls.

```bash
output_dir=~/nanomod_course_outputs/human
mkdir -p $output_dir
```

We can basecall and modification call our reads with 5mC methylation.
NOTE: the `-b 10 -c 1000` are internal `dorado` parameters which we've chosen to fit
our virtual machines.

```bash
input_dir=~/nanomod_course_data/human
model_files=~/nanomod_course_references/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v3.5.2
output_mod_bam=~/nanomod_course_outputs/human/PAM63974_pass_58881fec_60.twenty_random_reads.mod.bam
dorado basecaller $model_files $input_dir \
  --verbose -x cpu -b 10 -c 1000 --modified-bases 5mCG | \
    samtools view --threads 8 -O BAM -o $output_mod_bam
```

To call both 5hmC and 5mC methylation, use the flag `--modified-bases 5mCG_5hmCG`.
We are not going to call two modifications here.

We can view the output modification data in tabular form with `modkit`.

```bash
input_mod_bam=~/nanomod_course_outputs/human/PAM63974_pass_58881fec_60.twenty_random_reads.mod.bam
output_tsv=~/nanomod_course_outputs/human/PAM63974_pass_58881fec_60.twenty_random_reads.mod.bam.tsv
samtools index $input_mod_bam
modkit extract $input_mod_bam $output_tsv
```