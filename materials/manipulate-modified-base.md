---
layout: page
element: notes
title: Manipulation of base modification data
---

<!--
types of manipulation: thresholding, windowing, subsetting, pileup.
subset can be randomly, by genomic location, by mean density
modQC: how many modified, unmodified, no standard package so we'll have to do the tools ourselves.
three choices - our own program, samtools or modkit.
-->

<!-- TODO -->
```bash
./sample.test.awk sample.test.sam | samtools view -b | bedtools bamtobed -i - -tag XC
```
<!-- after making mod counts, can use commands like these  -->
```bash
samtools view -e '[XC]/qlen>0.02' -b -o testo.bam sample.bam
```

<!-- TODO: can introduce modbedtools https://github.com/lidaof/modbedtools 
https://doi.org/10.1016/j.xgen.2023.100455 -->
