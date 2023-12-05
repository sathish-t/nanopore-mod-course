---
layout: page
title: Computer Setup
---

## Overall Instructions

### If you are a participant in the training course

You will be given login details to virtual machines by the organisers or the trainers.
All the software and data required for the course is on these.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).

### If you are a self-study student

We use the following software packages in the course. 
We have tested the materials using the corresponding versions of the software.
NOTE: As the course material is under active development, we may not use all the software
below for the course.

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
```

#### Software packages used in python and R

We installed the following software packages on top of the base installs of python and R.

In Python, we installed:
- h5py
- modbampy
- numpy
- scikit-learn
- scipy
- pandas 
- pysam

In R, we installed the following software packages.
- ggplot2
- hexbin
- devtools
- ggthemes