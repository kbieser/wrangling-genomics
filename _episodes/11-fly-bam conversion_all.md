---
title: "FlyCURE - Variant Calling Workflow: Converting, Sorting, and Indexing bam files"
teaching: 0
exercises: 0
questions:
- "How do I convert sam to bam files to sort and index the alignment?"
objectives:
- "Understand the steps involved in variant calling."
- "Describe the types of data formats encountered during variant calling."
- "Use command line tools to perform variant calling."
keypoints:
- "Bioinformatic command line tools are collections of commands that can be used to carry out bioinformatic analyses."
- "To use most powerful bioinformatic tools, you'll need to use the command line."
- "There are many different file formats for storing genomics data. It's important to understand what type of information is contained in each file, and how it was derived."
---

# SAM/BAM format
The [SAM file](https://genome.sph.umich.edu/wiki/SAM),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not
have time to go into detail about the features of the SAM format, the paper by
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

~~~
#!/bin/bash

# what I do:

# step 1: sam to bam conversion
# step 2: sort bam by name (samtools sort -n)
# step 3: fixmate on bam
# step 4: samtools sort by coordinate
# step 5: mark dupes +remove (samtools markdup -r)
# step 6: indexes the markdup bam
# step 7: calculates the stats on that file


# run me in the folder with sam.gz (~/data/FlyCURE/results/bwa_out)
# if you want to save literal commands to a log file, run me like
# ../scripts/bam_factory.sh > mylogfile.txt
# therefore,
# intermediate folder is defined below as $inter
# clean folder is defined below as $clean

inter='../intermediate_bams'
mkdir -p $inter

clean='../clean_bams'
mkdir -p $clean

#this loop converts sam to bam
for i in *.sam.gz; do
  echo "converting $i to bam file"
  prefix=$(basename $i .sam.gz)
  echo samtools view -b -o ${inter}/${prefix}.bam $i
  samtools view -b -o ${inter}/${prefix}.bam $i &
done
wait

# sorting by name loop
for i in *.sam.gz; do
  echo "name sorting $i"
  prefix=$(basename $i .sam.gz)
  echo samtools sort -n -o ${inter}/${prefix}.nsort.bam ${inter}/${prefix}.bam &
  samtools sort -n -o ${inter}/${prefix}.nsort.bam ${inter}/${prefix}.bam &
done
wait

# fix mate loop
for i in *.sam.gz; do
  echo "fixmate $i"
  prefix=$(basename $i .sam.gz)
  echo samtools fixmate -r -m ${inter}/${prefix}.nsort.bam ${inter}/${prefix}.fixmate.bam &
  samtools fixmate -r -m ${inter}/${prefix}.nsort.bam ${inter}/${prefix}.fixmate.bam &
done
wait

# re-sort by coordinate
for i in *.sam.gz; do
  echo "coordinate sorting $i"
  prefix=$(basename $i .sam.gz)
  echo samtools sort -o ${inter}/${prefix}.csort.bam ${inter}/${prefix}.fixmate.bam &
  samtools sort -o ${inter}/${prefix}.csort.bam ${inter}/${prefix}.fixmate.bam &
done
wait

# markdup & remove pcr duplicates loop
for i in *.sam.gz; do
  echo "markdup removing dupes $i"
  prefix=$(basename $i .sam.gz)
  echo samtools markdup -r ${inter}/${prefix}.csort.bam ${clean}/${prefix}.markdup.bam &
  samtools markdup -r ${inter}/${prefix}.csort.bam ${clean}/${prefix}.markdup.bam &
done
wait

# index loop
for i in *.sam.gz; do
  echo "indexing $i"
  prefix=$(basename $i .sam.gz)
  echo samtools index ${clean}/${prefix}.markdup.bam &
  samtools index ${clean}/${prefix}.markdup.bam &
done
wait

# what I do:
# run flagstat on an indexed bam file and write that stat report to a log

# sorting by name loop
for i in *.sam.gz; do
  echo $index $bam to $log
  prefix=$(basename $i .sam.gz)
  echo samtools flagstat ${clean}/${prefix}.markdup.bam \> ${clean}/${prefix}.markdup.flags.log &
  samtools flagstat ${clean}/${prefix}.markdup.bam > ${clean}/${prefix}.markdup.flags.log &
done
~~~
{: .bash}

~~~
$ chmod +x bam_factory.sh
$ cd ~/data/FlyCURE/results/bwa_out
$ ../../scripts/bam_factory.sh > bam_factory.log
~~~
{: .bash}

We can watch it's progress by using the following command. Except this to run for a day.

~~~
tail -F bam_factory.log
~~~
{: .bash}
