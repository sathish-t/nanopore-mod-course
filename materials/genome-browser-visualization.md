---
layout: page
element: notes
title: Genome browser visualization of base modification data
---

In this session, we will visualize modification data in mod BAM files (see fig. below) using (1) genome
browsers (left) where we can rapidly scan data visually across different reads and different
regions on the genome, and (2) custom scripts (right) which allow us to see modification
density versus coordinate one read at a time.
We will also do some analysis on mod BAM files using `samtools` and `modkit` to show
how to re-construct the steps performed by genome browsers so that we can incorporate
it in our own scripts/commands.

![Compare IGV and rain plot visualization](compare_igv_and_rain.png)

Either visualization has its strengths and weaknesses, which is why it is good to know how to do both.
Genome browsers give you an overall idea and one can zoom in to regions on the genome
and see data across multiple reads. But, beyond recognizing which regions are highly
modified on a read (green) and which regions are not (grey), one cannot pick out any
details of modification density such as gradients as there is no interpolation
between the two colours to show any intermediary densities.
On the right, we have plotted a read using
a custom script, which shows raw modification data (grey) and windowed modification
data (red). Here, we can see details per read but we cannot see multiple reads at
the same time. We will explore the details of these visualizations and how to make
them in this session.

Modification pipelines have very similar steps till the modifications are called, and very different
steps after the modifications have been called.
Our goal in these sessions is to introduce several tools
that can be utilized to analyze modifications after they are called and show some sample ways
in which they can be executed.
The experimentalist must decide whether to use these tools, how to run these tools,
and/or whether new tools are needed depending on their experiment. 

## Visualizing modification calls with IGV

### Loading a mod BAM file

Let us open IGV on the virtual machines using the instructions from the figure below,
like we did in the [alignment]({{ site.baseurl }}/materials/sequence-align-pycoqc) session.
If you are a self-study student, please open IGV on your computer.

![Instructions on how to open IGV](open_igv.png)

Let us load the sacCer3 genome and the mod BAM file produced from the yeast dataset in the
[modification calling]({{ site.baseurl }}/materials/base-mod-detection) session following
the instructions in the figures below.

![Instructions on how to load genome](igv_load_genome_screenshot.png)
![Instructions on how to load BAM file](igv_load_file_screenshot.png)

You should see two tracks immediately below the reference genome on top.

### Colouring tracks by modification and exploring the visualization

Now, we have to choose the option of 'Color by modification'.
Please follow the instructions in the figure below.

![Instructions on how to colour by modification](igv_colour_by_mod.png)

Now, you can zoom in to the genome. Select any region of size around 10 to 50 kb.
You should see something similar to the figure below.

![IGV view with pileup annotation](igv_overall_view_with_pileup_annotated.png)

The gray and green tracks on top result from a pile-up analysis and
show the modified and unmodified count versus genomic position
(any averaging calculation that operates over all reads at every base is called a
pile-up analysis).

Try clicking on a read and see if you can get the following dialog boxes.

![IGV view with read dialog boxes](igv_get_read_details.png)

Try zooming in to a length scale of just a few bases. You should see individual
bases marked as modified and unmodified.

![IGV view with individual modified bases](igv_individual_modified_bases.png)

We can look at the reads for a few more minutes to get familiar with the visualization.
Pick a read id of a read that looks interesting to you and record it somewhere.
When we visualize single molecules, you can visualize this molecule.

## Reproducing IGV analyses with samtools

Before we proceed to single-molecule visualization, let us discuss how to do some of the analyses
above on the command line with `samtools` or `modkit`. This knowledge will come in handy when
you want to write scripts on your own.

### Subsetting mod BAM files with samtools

The first operation we will look at is subsetting the mod BAM file.
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
not just the part of the read corresponding to the region.
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