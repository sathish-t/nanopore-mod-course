---
layout: page
title: Computer Setup
---

## Instructions

### If you are a participant in the training course

You do not need to read any further than this section.
You will be given login details to virtual machines by the organisers or the trainers.
You will just need to download the latest version of the course repository
and create a few directories (instructions follow after this paragraph).
All other software and data required for the course is on the virtual machines, and
the commands to execute each piece of software will work straight away
on your Linux command line or on the GUI after you log in.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).

We need the latest version of the course repository.

```bash
cd ~/nanomod_course_scripts # go to the folder where scripts are stored
git clone  --depth 1 {{ site.github.repo }} # get latest version of the course repository
```

We store input data, scripts, references
and outputs in the directories `~/nanomod_course_data`, `~/nanomod_course_scripts`,
`~/nanomod_course_references` and `~/nanomod_course_outputs` respectively.
Please make them if they do not exist using the `mkdir` command.

### If you are a self-study student

We run our programs on the Linux command line and inspect the outputs either on the command
line or in the Linux desktop. We use the software packages below.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).
Please note that you are not expected to know how to code in python or R.
You would only need to be able to run programs written in these languages.

```bash
bedtools -version
# bedtools v2.29.2
DNAscent --version | head -n 1
# Version: 2.0.2
dorado --version
# 0.3.4+5f5cd02
guppy_basecaller --version | head -n 1
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 5.0.7+2332e8d
minimap2 --version
# 2.24-r1122
modbamtools --version
# modbamtools, version 0.4.8
modkit --version
# mod_kit 0.2.3
pod5 --version
# Pod5 version: 0.3.2
pycoQC --version
# pycoQC v2.5.2
python --version
# Python 3.10.2
samtools --version | head -n 2
# samtools 1.18
# Using htslib 1.18
R --version | head -n 1
# R version 4.1.2 (2021-11-01) -- "Bird Hippie"
aws --version
# aws-cli/2.14.5 Python/3.11.6 Linux/5.4.0-167-generic exe/x86_64.centos.7 prompt/off
git --version
# git version 2.25.1
wget --version | head -n 1
# GNU Wget 1.20.3 built on linux-gnu.
gunzip --version | head -n 1
# gunzip (gzip) 1.10
tar --version | head -n 1
# tar (GNU tar) 1.30
```

We need IGV v2.16.2 or later, which can be downloaded [here](https://igv.org/download/html/download.html).

We need the [DNAscentTools](https://github.com/DNAReplicationLab/DNAscentTools/) repository.
Get the latest version using the command below and then checkout to the version denoted by the commit
string `6bacc1f`.

```bash
mkdir -p ~/nanomod_course_scripts
# we store scripts in the above folder.
cd ~/nanomod_course_scripts
git clone https://github.com/DNAReplicationLab/DNAscentTools.git
cd DNAscentTools/
git checkout 6bacc1f
```

We need the latest version of the course repository.

```bash
mkdir -p ~/nanomod_course_scripts
cd ~/nanomod_course_scripts
# we store scripts in the above folder.
git clone --depth 1 {{ site.github.repo }}
```

#### Directory structure

We use the following four directories to store input data, scripts,
references, and outputs respectively.
You can use different directories if you wish.
Please substitute suitably throughout the course material
if you use different directories.

```bash
mkdir -p ~/nanomod_course_data
mkdir -p ~/nanomod_course_scripts
mkdir -p ~/nanomod_course_references
mkdir -p ~/nanomod_course_outputs
```

#### Software packages used in python and R

We need the following software packages on top of the base installs of python and R.

In Python, we need:
- h5py
- modbampy
- numpy
- scikit-learn
- scipy
- pandas 
- pysam

In R, we need:
- ggplot2
- hexbin
- devtools
- ggthemes


## Descriptions of software packages and links to their documentation and installation instructions

We list the software packages we use with a brief description and useful links below.

## Dorado  

Oxford NanoporeTech's latest, open-source basecaller.
Takes raw nanopore currents as input and produces DNA sequences and optionally modification calls as output
using a set of learned model parameters that associates bases with current characteristics.
Input formats are fast5 or pod5. Output format is .bam by default and other formats if requested.
- [Installation and documentation](https://github.com/nanoporetech/dorado)


## Modkit

Oxford NanoporeTech's powerful tool to convert modification data from a machine-readable format stored in BAM files
to a human-friendly, tab-separated format after optional, user-requested data processing steps. The tool can perform
sophisticated analysis like [differential methylation scoring](https://nanoporetech.github.io/modkit/intro_dmr.html).

- [Installation](https://github.com/nanoporetech/modkit)
- [Documentation](https://nanoporetech.github.io/modkit/) 