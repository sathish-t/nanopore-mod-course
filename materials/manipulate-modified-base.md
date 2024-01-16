---
layout: page
element: notes
title: Manipulation of base modification data
---

In this session, we will look at commands that manipulate modification 
data in mod BAM and produce summary statistics or plots.
Some of these commands are similar to those run behind the scenes
in a visualization software.
By running these commands ourselves, we get quantitative output
in a tabular format instead of an image format, and we can manipulate
this data further in our own scripts/pipelines.

The operations we will be looking at are:
- thresholding
- making histograms of modification probabilities
- subsetting (either randomly, or by genomic location, or by mean density)
- windowing modification calls
- pileup of modification calls 

![List of manipulations with mod BAM](manipulate_mod_bam.png)

We can achieve most of these operations using `modkit` and `samtools`.
For windowing and for measuring mean modification densities across reads,
we will use our own tools due to an absence of readymade, pre-existing tools.
In a later [session]({{ site.baseurl }}/materials/single-molecule-visualization),
we will use `modbamtools` for further manipulation.

## Refresher: Examine mod BAM files with `modkit extract`

As we have seen in a previous [session]({{ site.baseurl }}/materials/mod-bam-format),
the raw mod BAM format is pretty hard to parse for a person.
It is much more convenient to convert mod BAM files into a tabular format using
the `modkit extract` command.

```bash
input_bam_file= # fill with a suitable file
modkit extract $input_bam_file -
# - means the output goes to the standard output i.e. the screen.
# if your file is too large, pipe output through shuf and head.
# e.g. modkit extract $input_bam_file - | shuf | head -n 20
```

Please pay close attention to the first, second, third, fourth, and eleventh columns
which contain data corresponding to the read id, position along the read,
position along the reference, contig on the reference, and modification probability.

## Simple thresholding of mod BAM files with `modkit`

Thresholding is the process of converting soft modification calls to hard calls
i.e transforming a probability of modification per base per sequenced strand to a binary yes/no
(or a ternary yes/no/no etc. when multiple modifications are present).
In reality, a base on a DNA strand can be either modified or unmodified.
So, the probability of modification associated with every base in a mod BAM file
is due to uncertainty in the measurement procedure and for some applications,
this uncertainty is not important.
A thresholding step gets binary calls by converting calls above a probability
to 'yes' and calls below a probability to 'no', excluding bases where information
is missing.
The simplest threshold is 50% or 0.5.
We can execute the thresholding step using `modkit`.

```bash
input_mod_bam= # fill suitably
output_mod_bam= # fill suitably
modkit call-mods --no-filtering $input_mod_bam $output_mod_bam
```

Run `modkit extract` on the input and output mod bam files.
You will see that the modification probability is now either
a value very close to zero or a value very close to one.

The reason we do not get zero or one is because mod BAM probabilities
are discretized by 1/256. So, modkit outputs the midpoint of
the first and the last interval as the lowest and the highest
probabilities respectively.

### (optional) Multiple modifications

Please note that if there are multiple modifications, a base
can exist in more than two states: unmodified, modification of the first type,
modification of the second type etc. 
Now `call-mods` assigned the state with the highest measured probability to the base.

### Syntax for thresholds other than 50%

If we want to use a threshold other than 0.5, the syntax is 

```bash
input_bam= # fill suitably
output_bam= # fill suitably
mod_code= # fill with mod code. e.g.: T for our BrdU data, m for 5mC.
threshold= # fill with a number between 0 and 1.
modkit call-mods --mod-threshold $mod_code:$threshold --filter-percentile 0 $input_bam $output_bam
```

### (optional) Use two thresholds and move in-between bases to the missing category

In our procedure thus far, bases with probabilities above and below a threshold
are called as modified and unmodified respectively, and missing bases are left untouched.
Another approach uses two thresholds and calls bases above the higher threshold as modified,
calls bases below the lower threshold as unmodified,
and moves the bases between the two thresholds to the missing category.
The reasoning is that when a model assigns a modification probability close to
zero or 1, the model is confident in its call. When a model assigns a probability
close to 0.5, the model is not confident at all and is no better than the random baseline,
which the practitioners of this approach regard as no different from missing information.
This thresholding is accomplished through the `--filter-percentile` and/or
the `--filter-threshold` options of `modkit call-mods` in addition to the `--mod-threshold`
parameter. We will not be discussing this further; you can refer to the modkit documentation
[here](https://nanoporetech.github.io/modkit/advanced_usage.html#call-mods) if you are interested.


## Histograms of modification probabilities using `modkit`

In this section, we will calculate histograms of per-base modification probabilities.
One of its uses is helping us choose a threshold for a thresholding step that we were
discussing above.
Another use is to see if our experiment worked and we are seeing an expected number
of bases with modifications.

We will use the `modkit sample-probs` command. 
We operate the tool in the histogram mode `--hist` and want ten
bins in the histogram (`--buckets 10`), and wish to sample
1000 reads (`-n 1000`) to measure statistics (we are using 1000
to save runtime; usually a fraction of the BAM file is sufficient
for these calculations and we leave it to the experimenter to
determine what that fraction is for their dataset).

```bash
input_bam_file= # fill this suitably
modkit sample-probs $input_bam_file -o ./histogram --hist --buckets 10 -n 1000
```

The output files are deposited in the `./histogram` folder.
You should see three output files: `probabilities.txt`, `probabilities.tsv`
and `thresholds.tsv`.

<!-- TODO: fill what these files mean. what do we do about thresholds? -->

## Subsetting mod BAM files with samtools

Subsetting means we take a mod BAM file as input and produce a mod BAM file
as output by accepting or rejecting each read based on some criterion.

Whenever you zoom in to a region or click on a read in IGV, you are basically
subsetting a mod BAM file by region or by read id.
We will also discuss random subsetting which comes in handy when a BAM file
is very large and we are interested in a calculation whose result depends
very weakly on the number of reads, so it is sufficient to run the calculation
on a randomly-selected subset of reads.

Before we perform any subset, we first count the total number of reads we have
in a mod BAM file.

Count number of reads
```bash
input_file= # fill with whatever input file you want to use
samtools view -c $input_file
```

Subset by region:
Note that the subset will pick out entire reads that pass through a given region,
not just the part of the read corresponding to the region. Any subset <!-- TODO -->
```bash
contig=chrII
start=80000
end=90000
input_file= # fill with whatever input file you want to use
output_file= # fill with an output file name
samtools view -b -o $output_file $input_file $contig:$start-$end # perform subset
```

Let's do a quick check that our subset worked by counting the number of reads
and by examining their coordinates.

```bash
samtools view -c $input_file      # count reads of input file
samtools view -c $output_file     # count reads of output file
bedtools bamtobed -i $output_file | shuf | head -n 10 
    # look at a few output coordinates
    # you can run the bedtools command without the shuf and the head if
    # the output__file is only a few lines long.
```

Subset by read id:
```bash
# fill the following values. 
# use any suitable mod BAM file and some read id of interest you have recorded.
input_file=
read_id=
output_file=
samtools view -b -e 'qname=="'$read_id'"' -o $output_file $input_file # perform subset
samtools view -c $output_file # count reads
```

Get a random subset:
```bash
# fill the following values suitably. 
# use any suitable mod BAM file
input_file=
output_file=
fraction=0.05
samtools view -s $fraction -b -o $output_file $input_file # perform subset
samtools view -c $output_file # count reads
```

### Ways to perform pileup with samtools

<!-- TODO: flesh this out more -->

Get coverage
`samtools coverage dnascent.detect.mod.sorted.bam`

Get modification
`samtools mpileup -M dnascent.detect.mod.sorted.bam`

### Pileup with modkit

<!-- TODO: flesh this out more -->
<!-- TODO: do we need --no-filtering? -->
`modkit pileup --no-filtering --mod-thresholds T:0.03 dnascent.detect.mod.sorted.bam test.08jan24.bed`


<!--
introduce thresholding and windowing but need not get into details here.
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

<!-- TODO: 
modkit sample-probs -n 1000 --out-dir ./test --hist --buckets 20 gm12878_ul_sup_megalodon_HP_chr20.bam
-->