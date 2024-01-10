#! /usr/bin/awk -f
# Written using ChatGPT.
# Usage: samtools view -h sample.bam | awk -f count_mods_per_read.awk
# Input: Pass thresholded mod BAM files with only one type of modification as input
# Output: Plain text, same as input mod BAM but with new XC tag per line
#         with modification count of that line
BEGIN{OFS=FS="\t"}
!/^@/ {
  found=0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^ML:B:C,/) {
      split($i, arr, ",");
      count=0;
      for (j=2; j<=length(arr); j++) {
        if (arr[j] == 255) count++
      }
      $NF = $NF "\tXC:i:" count; found=1
    }
  }
  {
    if (!found) $NF = $NF "\tXC:i:0"
  }
  found=0
} 1