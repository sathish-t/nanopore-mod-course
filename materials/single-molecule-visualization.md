---
layout: page
element: notes
title: Visualisation of base modification data - part II
---

## Pileup of reference-anchored mod BAM files with `modkit`

A pileup is any calculation that produces one number per base on a reference
genome by performing an operation across all data available at that base
across all reads passing through that base.
In addition to the percent-modification track produced by modbamtools
above, we have also encountered pileups in the previous
[session]({{ site.baseurl }}/materials/genome-browser-visualization)
on visualization in genome browsers as shown in the figure below.

![IGV view with pileup annotation](igv_overall_view_with_pileup_annotated.png)

Per base on the reference genome, we can calculate
- the total number of reads
- the total number of modifications
- the fraction of modified reads

```bash
input_mod_bam=         # fill suitably
output_dir=            # fill suitably
modkit pileup --no-filtering --mod-thresholds T:0.5\
  $input_mod_bam "$output_dir"/pileup.bed
```

Now, one can inspect a few lines from the output file

```bash
# inspect the first few lines
head -n 20 "$output_dir"/pileup.bed
# inspect a few randomly chosen lines
cat "$output_dir"/pileup.bed | shuf | head -n 20
```

### (optional) Coverage using `bedtools`

An alternate way to get the coverage,
which is a count of the number of reads passing
through each base, one can do

```bash
input_mod_bam=         # fill suitably
output_dir=            # fill suitably
mkdir -p "$output_dir" # make output directory if need be
bedtools genomecov -ibam $input_mod_bam -bga >\
  "$output_dir"/coverage.bedgraph
```

Now, one can inspect a few lines from the output bedgraph

```bash
# inspect the first few lines
head -n 20 "$output_dir"/coverage.bedgraph
# inspect a few randomly chosen lines
cat "$output_dir"/coverage.bedgraph | shuf | head -n 20
```

This method will work even if the BAM file has no modification
information as coverage is just a count of number of reads
passing through each position on the reference irrespective
of whether they contain modifications or not.

### (optional) Modification pileup with `samtools`

One can also perform pileups of modification with `samtools`.
The command is specified below.
The output format is a little hard to understand and we will not be discussing this further.
Please consult the documentation [here](https://www.htslib.org/doc/samtools-mpileup.html)
if you want to learn more.

```bash
input_mod_bam= # fill suitably
output_file= # fill suitably
samtools mpileup -M $input_mod_bam > $output_file
```