---
layout: page
element: notes
title: MOD BAM format
---

In this page, we are going to discuss how single-molecule modification data
is stored in BAM files using optional tags.
We are going to be building upon our knowledge of the BAM file format, so please
make sure you are familiar with it as explained in [this]({{ site.baseurl }}/materials/sequence-align-pycoqc)
session.

<!-- TODO:
 mention need not know this -->

We will illustrate multiple ways to represent modification information
using an example molecule with the sequence `TCGCCTAGCG`.
Let us assume the read has not been aligned to a reference genome for simplicity.

## Scenario 1: Bases are either modified or unmodified with no ambiguity or missing information

In the simplest scenario, we just want to mark which bases are modified and which are not per molecule.

Let us consider the example where the second and the third cytosines of the sequence of interest above have been
modified to 5-methylcytosine (5mC) and the others are unmodified.
In other words, when moving along cytosines on the sequence,
we have to skip `1` cytosine to get to the first modified base and then skip `0` cytosines to get
to the second modified base.
So, the modification tag reads `MM:Z:C+m,1,0;`.
In addition to the information that we need to execute 1 skip and 0 skips (`,1,0`), the tag tells
us that we are looking at cytosine modifications (`C`), the modification under
question is 5-methylcytosine (`m`), and that modification information is on the
same strand as basecalling (`+`). Put together, the BAM line might look like (spaces have replaced tabs for simplicity):

```text
1578e830-9886-46c8-977b-f0645dac040b 4 * 0 255 * * 0 0 TCGCCTAGCG * MM:Z:C+m,1,0;
```

Henceforth, we will drop the other fields in the BAM line and focus only on
the modification tags for simplicity. The first five characters of the tag
are always `MM:Z:`; older mod bam files may also use `Mm:Z:`.

### Standard modifications have one letter codes or numeric codes

We have used the one-letter code `m` to represent 5mC in the example above.
These methylation codes are a part of the standard BAM format
specification [here](https://samtools.github.io/hts-specs/SAMtags.pdf).
5-Hydroxymethylcytosine (5hmC) is represented with an `h` or the numeric
code `76792`.
So if our example above dealt with 5hmC modification, the tag would
read as `MM:Z:C+h,1,0;` or `MM:Z:C+76792,1,0;`.

If this sample molecule were a part of our yeast dataset where thymidines are substituted by
BrdUs (numeric code 472252), the tag would read `MM:Z:T+472552,0;` if we wanted
to record that the first thymidine in our sample sequence `TCGCCTAGCG` was substituted.
As older modification tools did not recognize numeric codes,
we have represented BrdU substitution using the ambiguous one letter code T
instead i.e. `MM:Z:T+T,0;`.

### (optional) Multiple types of modifications are present

Let us consider the scenario where both 5mC and 5hmC modifications are present in the sample sequence `TCGCCTAGCG`.
The tag might read (a space has replaced a tab for simplicity) `MM:Z:C+m,1,0; MM:Z:C+h,3;`.
The interpretation is that the second and third cytosines (execute one skip and zero skips) have been modified to 5mC and
the fourth (execute three skips) has been modified to 5hmC.

### (optional) Modifications are present on the basecalled strand and its complement

ONT have released a 'duplex' method where both a strand and its complement from the same cell are sequenced.
So, one can obtain modification status on both a strand and its complement.
In our sample sequence `TCGCCTAGCG`, let's say the second and third cytosines have been modified to 5mC,
and the first cytosine on the complementary strand which pairs with the highlighted G `TC_G_CCTAGCG` has
been modified to 5mC.
The tag would then read (replace spaces with tabs) `MM:Z:C+m,1,0; MM:Z:G-m,0;`.

## Scenario 2: Base modifications are specified using probabilities