---
layout: exercise
topic: BAM files
title: Exploring an unknown BAM file
---

We perform some explorations on the BAM file located at
`~/nanomod_course_data/human/bonito_calls.subset.sorted.bam`.

Tip: Any program that creates the BAM file or alters it in some way will
usually add a line to the BAM file header.
You can use `samtools view -H $mod_bam` (substitute `$mod_bam` suitably) to view the header.
If the header is many lines, you can pipe the output to `head` or `tail` or `shuf` to view
a few lines e.g. `samtools view -H $mod_bam | head -n 20`.

1. The modification calling was performed using a different program which is
very similar to `dorado`. Can you find the program and the exact command?

2. The course trainers subset a 'Cliveome' BAM file to make this file.
Can you locate the subset command and infer how the subset was performed?
This should tell you which region of the genome the BAM file maps to.

3. How many reads are in the BAM file?

4. How many secondary/supplementary reads are in the BAM file?