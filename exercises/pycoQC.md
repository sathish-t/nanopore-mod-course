---
layout: exercise
topic: pycoQC
title: Running pycoQC on BAM files
---

Run pycoQC on the BAM file `~/nanomod_course_data/yeast/subset_2.sorted.bam`
and the corresponding sequencing summary file
`~/nanomod_course_data/yeast/sequencing_summary.subset_2.txt`.
Do you see anything unusual about the coverage overview and the read length distribution?

```bash
# command syntax below.
# please fill files suitably
pycoQC -f $input_seq_sum -a $input_bam -o $output_html \
  -j $output_json --quiet
```