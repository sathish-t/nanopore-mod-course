#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J lPlotRead
#SBATCH --mail-type=END,FAIL
#SBATCH --time 0:19:59
#SBATCH --constraint="intel"

# preamble
# --------

# a sample program execution line and its meaning follows

# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir
# can use sbatch in place of bash.
# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir mash
# can also optionally provide a prefix like 'mash' in the example above
# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir mash annotation_file
# can also optionally provide a file called annotation_file in the line above.
# caution: if annotation_file is provided, then a prefix must be provided.
# The annotation file is space-separated and has three columns (without headers):
# start, end, label. start and end refer to locations on the reference genome,
# and label can be 'origin', 'termination', 'leftFork', 'rightFork', or 'pause'.

# this script is a limited version of plot_read.sh.
# the line above extracts data from the read with read id = readID from sample.mod.bam.
# Data is windowed in 300 thymidines after thresholding
# i.e. calling T as BrdU if probability > 0.5 (threshold cannot be changed in this script).
# If you do not want to show a windowed curve, pass a window size of 0 in the command line invocation.
# The plot and plot data are sent to output_dir and have names like
# plot_readID.png, plot_data_readID.
# If a suffix is specified, say 'mash', then the filenames are mash_readID.png, mash_data_readID
# Plot has two components: the raw data and the windowed data.

# if you want any other functionality, refer to plot_read.sh.

# stop execution if any command fails
set -e

# set output directory, making it if it doesn't exist
mkdir -p "$4"
op_dir=$(cd "$4"; pwd)

# load packages
pwd=$(pwd)
config_dir=..
cd "$config_dir"
source load_package.sh -R -python -samtools -bedtools

{
  # get information about the read
  bash get_information_from_read_id.sh "$1" "$2";

  # get raw data
  samtools view -b -e 'qname=="'"$2"'"' "$1" |\
    bedtools bamtobed -i stdin |\
    awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
    sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
    python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$1"

} > "$op_dir"/"${5:-plot}"_temp_"$2"

{

    # output raw data with associated windows of 1 base each
    < "$op_dir"/"${5:-plot}"_temp_"$2" awk '/^#/ {print} !/^#/ {print $1 " " $2 " " $2+1 " " $3 " rawDetect"}'

    # window data
    if ! [ "$3" -eq 0 ];
    then
      < "$op_dir"/"${5:-plot}"_temp_"$2" \
        sed '1i\detectIndex\tposOnRef\tval' |\
        python get_mean_brdU_window.py --window "$3" --thres 0.5 |\
        awk '{print $1 " " $3 " " $4 " " $2 " winDetect"}';
    fi

} > "$op_dir"/"${5:-plot}"_data_"$2"

# plot data
cd "$pwd"
if [ "${6:-flyingTurtleMoons}" == "flyingTurtleMoons" ];
then
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$op_dir"/"${5:-plot}"_data_"$2" "$op_dir"/"${5:-plot}"_"$2".png;
else
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$op_dir"/"${5:-plot}"_data_"$2" "$op_dir"/"${5:-plot}"_"$2".png\
    "$6";
fi

# delete temporary file
rm "$op_dir"/"${5:-plot}"_temp_"$2"