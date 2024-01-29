---
layout: exercise
topic: Modification analysis
title: Compare modification levels
---

In the commands we executed on the first day, we used
`DNAscent forkSense` to call features such as initiation
and termination sites on single molecules.
In the case study today, we learned that in the experimental
protocol used in the yeast dataset, sites that are replicated
early are more modified than sites that are replicated late.
In this exercise, please demonstrate that sequences
on single molecules corresponding to initiation sites
are more modified than sequences on single molecules
that correspond to termination sites.

Your calculation should use the following files as input:
- `~/nanomod_course_outputs/yeast/dnascent.detect.mod.sorted.bam`
- `~/nanomod_course_outputs/yeast/origins_DNAscent_forkSense.bed`
- `~/nanomod_course_outputs/yeast/terminations_DNAscent_forkSense.bed`

You can answer the question in any manner you like.
You can produce a quantitative answer or a qualitative answer.
We suggest that you formulate a plan and try to execute it.
If it turns out that your plan requires a lot of work, then
please abandon it and try to think up a simpler solution.

For some answers, you may want to use the bed files above as an
input to programs such as `modkit`.
Please be advised that although `DNAscent forkSense`
outputs files ending in `.bed`, these are not in a correct format.
So, if you intend to use these files as inputs, please convert
them to a bed format before doing so.
The following commands may be helpful for the conversion.

```bash
input_forksense_bed= # fill suitably
output_forksense_bed= # fill suitably
awk -v OFS="\t" '{print $1, $2, $3, $4, 100, $8=="fwd"?"+":"-"}' \
  $input_forksense_bed > $output_forksense_bed
```