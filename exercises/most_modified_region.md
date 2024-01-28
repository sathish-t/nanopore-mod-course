---
layout: exercise
topic: Modification analysis
title: Most modified region
---

The file `~/nanomod_course_data/yeast/subset_2.sorted.bam`
contains reads from the yeast dataset longer than 30 kb
that pass through the region chrVI:150000-250000.
Our goal is to identify which of the following three regions
is the most modified:
- chrVI:220000-230000
- chrVI:150000-160000
- chrVI:190000-200000 

You can solve this in whichever way you want and produce
a qualitative or a quantitative solution.
For the official solution, we want you to implement
the steps below:

- Our molecules fall into two categories: newly-synthesized
DNA with BrdU in it and DNA that was synthesized before the
current cell cycle and does not have BrdU.
In other words, some DNA strands contain
modifications and some do not. The ones that do not contain
modifications will still have some bases called as modified
due to the unavoidable false-positive rate of modification callers.
Form a subset of reads that contain only reads with modifications using
any reasonable procedure.
- Subsample the subset formed above if you desire so.
- Calculate modification amounts for each region
and identify the most modified region.