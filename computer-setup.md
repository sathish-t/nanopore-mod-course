---
layout: page
title: Computer Setup
---

## Instructions

### If you are a participant in the training course

You will be given login details to virtual machines by the organisers or the trainers.
All the software and data required for the course is on these.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).

### If you are a self-study student

We use the software packages below. Please install the corresponding versions.
You are not expected to know how to code in python or R.
You would only need to be able to run programs written in these languages.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).

```bash
bedtools -version
# bedtools v2.29.2
DNAscent --version | head -n 1
# Version: 2.0.2
dorado --version
# 0.3.4+5f5cd02
minimap2 --version
# 2.24-r1122
modbamtools --version
# modbamtools, version 0.4.8
modkit --version
# mod_kit 0.2.3
pycoQC --version
# pycoQC v2.5.2
python --version
# Python 3.10.2
samtools --version | head -n 2
# samtools 1.18
# Using htslib 1.18
R --version | head -n 1
# R version 4.1.2 (2021-11-01) -- "Bird Hippie"
git --version
# git version 2.25.1
```

We also need the [DNAscentTools](https://github.com/DNAReplicationLab/DNAscentTools/) repository.
Get the latest version using the command below and then checkout to the version denoted by the commit
string `6bacc1f`.

```bash
git clone https://github.com/DNAReplicationLab/DNAscentTools.git
cd DNAscentTools/
git checkout 6bacc1f
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