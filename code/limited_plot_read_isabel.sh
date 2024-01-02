#!/bin/bash

# Goal: to extract raw and windowed data from one read and to plot it.

# Usage: > bash limited_plot_read.sh <sample.mod.bam> <readID> <mod_code> <position> <threshold> <window_size> <output_dir>

# Script: limited_plot_read.sh

# Input file: sample.mod.bam

# ReadID: readID
#  id from the read of interest
#  e.g., 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575

# Modification code: mod_code
#  DNA base modification of interest
#  e.g., T for any T modification, h for 5hmC, m for 5mC, a for 6mA

# Position of DNA base modification: position
#  must be either not_use_ref or use_ref.
#  not_use_ref: the position is the value of the column forward_read_position
#  use_ref: the position is the value of the column ref_position

# Threshold for modification calling: threshold
#  e.g., a thres = 0.5 means that a base with a modification probability >= 0.5 is modified and < 0.5 is unmodified

# Window size: window_size
#  number of bases that are used to window the data after thresholding
#  e.g., a win = 300 means that data is windowed in 300 bases after thresholding

# Output directory: output_dir

# Output files:
#  plot_png: a  plot with the DNA base modification raw data and windowed data.
#  one_read_bam: a mod.bam file with the read of interest.
#  one_read_bam.tsv: a tsv file with the modification information from the read of interest.

# Check the number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: bash limited_plot_read.sh <sample.mod.bam> <readID> <mod_code> <position> <threshold> <window_size> <output_dir>"
    exit 1
fi

# Stop execution if any command fails
set -euxo pipefail

# Set output directory, making it if it doesn't exist
mkdir -p "$7"
op_dir=$(cd "$7"; pwd)

# Make a bam file with just the read of interest
one_read_bam="$op_dir"/"plot"_"$2".bam
samtools view -b -o $one_read_bam -e 'qname=="'"$2"'"' "$1"
samtools index "$op_dir"/"plot"_"$2".bam

# Extract modification information from the bam file
# NOTE: the --force overwrites any pre-existing output files.
modkit extract --force "$one_read_bam" "$one_read_bam".tsv

# Use python scripts to extract raw data and window it
mod_code="$3"
position="$4"
threshold="$5"
window_size="$6"

python extract_raw_mod_data.py $mod_code $position "$one_read_bam".tsv "$one_read_bam"_rawDetect.tsv
python window_mod_data.py $threshold $window_size "$one_read_bam"_rawDetect.tsv "$one_read_bam"_winDetect.tsv
cat "$one_read_bam"_rawDetect.tsv" > "$one_read_bam".data && tail -n +2 $one_read_bam"_winDetect.tsv >> "$one_read_bam".data 

# plot data
Rscript ./plot_one_read_w_win.R "$one_read_bam".data "$one_read_bam".png;
