#!/bin/bash

# explanation
# --------
# a sample program execution line and its meaning follows
# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir

# the line above extracts data from the read with read id = readID from sample.mod.bam.
# Data is windowed in 300 bases after thresholding
# i.e. calling a base as modified if probability > 0.5 (threshold can be changed in the script).
# The plot and plot data are sent to output_dir and have names like
# plot_readID.png, plot_data_readID.
# Plot has two components: the raw data and the windowed data.

# stop execution if any command fails
set -euxo pipefail 

# set output directory, making it if it doesn't exist
mkdir -p "$4"
op_dir=$(cd "$4"; pwd)

# make a bam file with just the read of interest
one_read_bam="$op_dir"/"plot"_"$2".bam
samtools view -b -o $one_read_bam -e 'qname=="'"$2"'"' "$1" 
samtools index "$op_dir"/"plot"_"$2".bam

# extract modification information from this file
# NOTE: the --force overwrites any pre-existing output files.
modkit extract --force "$one_read_bam" "$one_read_bam".tsv

# use python to extract raw data and window it
{
  python extract_raw_mod_data.py T "$one_read_bam".tsv 
  python window_mod_data.py T 300 "$one_read_bam".tsv
} > "$one_read_bam".data

# plot data
Rscript ./plot_one_read_w_win_or_model_if_needed.R "$one_read_bam".data "$one_read_bam".png;