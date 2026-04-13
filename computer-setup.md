---
layout: page
title: Computer Setup
---

## Instructions

### Style of computer commands used in the course

In the course material, we follow the style of naming all the input parameters
to each command. For example, let's say we want to see the contents of a file
with DNA data located at `~/somewhere/someplace/somedata.txt` using the `cat`
command. The usual style of running this command is

```bash
# do not run this command! It's for illustrative purposes only.
cat ~/somewhere/someplace/somedata.txt
```

However, in the course, we declare a variable called `dna_data`
first which contains the path to the file. Then, we pass this
variable as an input to the command as shown below.

```bash
# do not run this command! It's for illustrative purposes only.
dna_data=~/somewhere/someplace/somedata.txt
cat $dna_data
```

It is easier to see what the command is doing in our style as we have
a name associated with the file i.e. instead of thinking 'we are viewing
the contents of a file located at this path', we think 'we are viewing dna data'.

This becomes even more relevant
when we run commands with many input parameters; see the two styles of running
the same command below.

```bash
# style 1
# do not run this command! It's for illustrative purposes only.
pycoQC -f ~/nanomod_course_outputs/yeast/sequencing_summary.txt \
   -a ~/nanomod_course_outputs/yeast/aligned_reads.sorted.bam \
   -o ~/nanomod_course_outputs/yeast/pycoQC/analysis.html \
  -j ~/nanomod_course_outputs/yeast/pycoQC/analysis.json --quiet
```

```bash
# style 2
# do not run this command! It's for illustrative purposes only.
input_seq_sum=~/nanomod_course_outputs/yeast/sequencing_summary.txt
input_bam=~/nanomod_course_outputs/yeast/aligned_reads.sorted.bam
output_html=~/nanomod_course_outputs/yeast/pycoQC/analysis.html
output_json=~/nanomod_course_outputs/yeast/pycoQC/analysis.json
pycoQC -f $input_seq_sum -a $input_bam -o $output_html \
  -j $output_json --quiet
```

This is a question of style. Although we prefer style 2 above,
you can execute the command in whichever way you want.

### If you are a participant in the Earlham Institute training course

You do not need to read any further than this section.
You will be given login details to virtual machines by the organisers or the trainers.
The trainers will show you how to login and how to access a terminal and a graphical
user interface (GUI) on the virtual machines, and how to paste commands if you want
to copy and paste commands from somewhere.
All software and data required for the course is on the virtual machines, and
the commands to execute each piece of software will work straight away
on your Linux command line or on the GUI after you log in.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing-2025).

We store input data, scripts, references
and outputs in the directories `~/nanomod_course_data`, `~/nanomod_course_scripts`,
`~/nanomod_course_references` and `~/nanomod_course_outputs` respectively.

#### Following along in the course

In the hands-on sessions, you should see a screen with a side-by-side view of
this course website and a terminal or a GUI window in the
video from the videoconferencing software.
The instructor will execute commands from this course website
in the terminal or the GUI window.
You can log in to a virtual machine and type commands along with the instructor.

Some of our commands are long and go over several lines.
Please make sure to maximize your terminal window to make these commands easy to view.

In the course material, we do not always write 'Run this command' before
every command block. Please follow what the speaker is doing or make a judgement call
about whether a block needs to be run.

### If you are a self-study student

We run our programs on the Linux command line and inspect the outputs either on the command
line or in the Linux desktop. We use the software packages below.
Please also consult the course prerequisites listed [here](https://www.earlham.ac.uk/events/detection-dna-base-modification-using-nanopore-sequencing).

This repository contains two files to help you get set up:
- [`software.def`]({{ site.github.repo }}/blob/main/software.def)
  — a Singularity definition file that builds a container image
  with all the software packages listed below
  (except DNAscent, which has its own
  [Singularity image](https://cloud.sylabs.io/library/mboemo/dnascent/dnascent)).
- [`.bashrc`]({{ site.github.repo }}/blob/main/.bashrc)
  — a shell configuration file that defines convenient wrapper
  functions so that each command (e.g. `samtools`, `modkit`)
  transparently runs inside the Singularity container.

To build the main software image and pull the DNAscent image,
run the following commands.
Please note that the software packages may require around 10 GB
of disk space.

```bash
sudo singularity build software.img software.def
singularity pull DNAscent.sif library://mboemo/dnascent/dnascent:4.1.1
```

```bash
bedtools -version
# bedtools v2.31.1
DNAscent --version | head -n 1
# Version: 4.1.1
dorado --version
# 0.9.6+0949eb8d
minimap2 --version
# 2.30-r1287
modbamtools --version
# modbamtools, version 0.4.8
modkit --version
# modkit 0.6.1
nanalogue --version
# nanalogue 0.1.9
pod5 --version
# Pod5 version: 0.3.36
pycoQC --version
# pycoQC v2.5.2
samtools --version | head -n 2
# samtools 1.23.1
# Using htslib 1.23.1
aws --version
# aws-cli/2.34.15
jq --version
# jq-1.6
# The following are system tools. Any reasonable version will suffice.
# They run from the host OS, not the container.
git --version
wget --version | head -n 1
gunzip --version | head -n 1
tar --version | head -n 1
```

We need IGV v2.16.2 or later, which can be downloaded [here](https://igv.org/download/html/download.html).

We need nanalogue-gui v0.2.7 or later, which can be downloaded [here](https://github.com/sathish-t/nanalogue-gui/releases).

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

## Descriptions of software packages and links to their documentation and installation instructions

We list the software packages we use with a brief description and useful links below.
Please install the versions we have listed [above](#if-you-are-a-self-study-student).

## Dorado  

Oxford NanoporeTech's latest, open-source basecaller.
Takes raw nanopore currents as input and produces DNA sequences and optionally modification calls as output
using a set of learned model parameters that associates bases with current characteristics.
Input formats are fast5 or pod5. Output format is .bam by default and other formats if requested.
- [Installation and documentation](https://github.com/nanoporetech/dorado)

## Nanalogue

Tool to parse and analyse BAM/Mod BAM files with a single-molecule focus.
Extracts and processes information from BAM files, with a particular focus on
single-molecule aspects and DNA/RNA modifications.
- [Installation and documentation](https://github.com/DNAReplicationLab/nanalogue)
- [Cookbook](https://www.nanalogue.com)

## Nanalogue-gui

Electron GUI for interactive sequence data analysis and curation with a focus
on single molecules and DNA/RNA modifications.
Provides QC, Swipe (annotation curation), Locate Reads, and AI Chat modes.
- [Installation and documentation](https://github.com/sathish-t/nanalogue-gui)

## Modkit

Oxford NanoporeTech's powerful tool to convert modification data from a machine-readable format stored in BAM files
to a human-friendly, tab-separated format after optional, user-requested data processing steps. The tool can perform
sophisticated analysis like [differential methylation scoring](https://nanoporetech.github.io/modkit/intro_dmr.html).

- [Installation](https://github.com/nanoporetech/modkit)
- [Documentation](https://nanoporetech.github.io/modkit/) 

## DNAscent

Reference-anchored modification-calling program. Uses raw nanopore currents, the reference genome, and
alignment information to output a probability of BrdU substitution at each thymidine in a sequenced strand.
Uses a set of learned model parameters that help distinguish between currents produced by thymidines and
BrdUs in different k-mer environments. 
- [Installation](https://github.com/MBoemo/DNAscent)
- [Documentation](https://dnascent.readthedocs.io/en/latest/base.html)

## Minimap2

Alignment software that finds the best-fit location of a given sequence on a linear reference genome.
- [Installation and documentation](https://github.com/lh3/minimap2)

## pod5

ONT-written package for reading and writing files in the pod5 format.
- [Installation and documentation](https://github.com/nanoporetech/pod5-file-format)

## jq

jq is a software used on the command line to process JSON data.
- [Homepage](https://jqlang.org)

## Other software

We will discuss the other software packages as and when we encounter them in the course.
- [Installation and documentation for modbamtools](https://rrazaghi.github.io/modbamtools/)
- [Installation and documentation for samtools](https://www.htslib.org)
- [Installation and documentation for bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Installation and documentation for pycoQC](https://github.com/a-slide/pycoQC)
