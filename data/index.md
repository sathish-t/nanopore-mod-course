---
layout: page
title: Datasets
---

## Instructions

We will be working with two datasets in this course: a yeast DNA dataset
with some thymidines substituted by BrdU and a human dataset of
methylated DNA.
The yeast dataset is at this [link](https://ckan.earlham.ac.uk/dataset/nat-meth-2019-subset-for-nanomod-course)
and the human dataset is at this [link](https://labs.epi2me.io/cliveome_5mc_cfdna_celldna/).
Please read on to see how to download the datasets and prepare them for use.

### If you are a participant in the Earlham Institute training course

You do not need to read any further than this section.
You will be given login details to virtual machines by the organisers or the trainers.
The input data you need has already been downloaded onto your computer.
You need to copy it into the correct location.
Please execute the following commands

```bash
mkdir -p ~/nanomod_course_data
cd ~/nanomod_course_data
cp -r /mnt/trainings/BaseMod2024/* .
```

You should now see two directories called `yeast` and `human`.
Data preparation is now complete.

### If you are a self-study student

#### Yeast

Please download the yeast dataset from this
[link](https://ckan.earlham.ac.uk/dataset/nat-meth-2019-subset-for-nanomod-course).
Please put it in a suitable folder and untar it using the command 

```bash
filename= # substitute filename suitably
tar -xzvf $filename
```

Then, move the contents of the resultant `for_ckan` folder to a suitable location.
In the course, we place it under `~/nanomod_course_data/yeast`, but you
can use any suitable location.
Just remember to use the correct paths when executing the commands in each session.

#### Human

The human dataset we are using is > 1 TB but we will be using only a small fraction of it.
As any human dataset is of a sensitive nature, we are unable to provide you the subset
we are using directly and you will have to recreate it using the commands below.
You will need to download tens of GB and use tens of minutes of computational time.

Please make a directory called `~/nanomod_course_data/human` to store the data.
As stated in the yeast subsection above, you can use any other directory.
Just remember to use the correct paths when executing commands in the course.

We first make a subset of a mod BAM file as shown below.
```bash
cd ~/nanomod_course_data/human
mod_bam=http://ont-open-data.s3.amazonaws.com/cliveome_kit14_2022.05/gdna/basecalls/PAM63974/bonito_calls.bam
samtools view -h -b -e 'rlen>=30000'  $mod_bam chr20:58000000-60000000  > bonito_calls.subset.bam
```

Using `samtools sort` and `samtools index`,
sort and index the file to form `bonito_calls.subset.sorted.bam` and `bonito_calls.subset.sorted.bam.bai`.
You will learn how to run these commands on day 1, so if you do not know how to run them,
just come back to this section after day 1.

Next, we make a random subset of nanopore currents from a fast5 file after converting it to pod5.

```bash
cd ~/nanomod_course_data/human
wget http://ont-open-data.s3.amazonaws.com/cliveome_kit14_2022.05/gdna/flowcells/ONLA29134/20220510_1127_5H_PAM63974_a5e7a202/fast5_pass/PAM63974_pass_58881fec_60.fast5
pod5 convert fast5 ./PAM63974_pass_58881fec_60.fast5 --output ./PAM63974_pass_58881fec_60.pod5
pod5 view ./PAM63974_pass_58881fec_60.pod5 --ids --no-header | shuf | head -n 20 > twenty_read_ids.txt
pod5 filter ./PAM63974_pass_58881fec_60.pod5 --ids twenty_read_ids.txt --output ./PAM63974_pass_58881fec_60.twenty_random_reads.pod5
```

Next, we run pycoQC to judge the quality of the dataset.
In the course run by the Earlham Institute, we run the program ourselves and
give the results to the course participants to minimize computer runtime.
Unfortunately, we cannot host these files on the internet and give them to you.
The commands below require downloading tens of GB of data
and tens of minutes of computational time.

```bash
cd ~/nanomod_course_data/human
wget http://ont-open-data.s3.amazonaws.com/cliveome_kit14_2022.05/gdna/flowcells/ONLA29134/20220510_1127_5H_PAM63974_a5e7a202/sequencing_summary_PAM63974_58881fec.txt
wget http://ont-open-data.s3.amazonaws.com/cliveome_kit14_2022.05/gdna/basecalls/PAM63974/bonito_calls.bam 
wget http://ont-open-data.s3.amazonaws.com/cliveome_kit14_2022.05/gdna/basecalls/PAM63974/bonito_calls.bam.bai 
pycoQC -f sequencing_summary_PAM63974_58881fec.txt\
  -a bonito_calls.bam -o ./analysis.html -j ./analysis.json
```

Data preparation is now complete.