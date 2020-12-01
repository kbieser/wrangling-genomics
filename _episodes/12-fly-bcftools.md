---
title: "FlyCURE - Variant Calling Workflow: bcftools"
teaching: 0
exercises: 0
questions:
- "How do I identify unique SNPs"
objectives:
- "Call unique SNPs by genomic coordinate for each mutant."
- "Write for loops for bcftools mpileup and calls."
keypoints:
- "Bioinformatic command line tools are collections of commands that can be used to carry out bioinformatic analyses."
- "To use most powerful bioinformatic tools, you'll need to use the command line."
- "There are many different file formats for storing genomics data. It's important to understand what type of information is contained in each file, and how it was derived."
---

**Launch the app from username/data and RUN CELL ONE OF YOUR PERSISTANCE NOTEBOOK before starting.**

# Review the markdup files

Let's first review the outputs from our bam_factory script.

~~~
$ cd ~/data/FlyCURE/results/clean_bams
$ less A44.markdup.flags.log
~~~
{: .bash}

~~~
34957207 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
33092 + 0 supplementary
0 + 0 duplicates
34957207 + 0 mapped (100.00% : N/A)
34766218 + 0 paired in sequencing
17383109 + 0 read1
17383109 + 0 read2
34250592 + 0 properly paired (98.52% : N/A)
34766218 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
364870 + 0 with mate mapped to a different chr
123478 + 0 with mate mapped to a different chr (mapQ>=5)
~~~
{: .output}

This output is the final alignment statistics after completing all of the prior samtools scripts. From the first line, you can see that ~35 million reads passed quality control and 0 reads failed. If you continue down the list, 100% of the reads mapped to the genomic reference and 98.52% of both read 1 and read 2 paired and mapped properly. This tells us we can have high confidence in our alignment. If you review each of the .markdup.flags.log files, you should see that each of them have equally high percentages.

# Run bcftools for making base calls

We now have the full alignments after running samtools. The next step is to use a program called bcftools.  Bcftools will take the alignment files and produce a list of unique SNPs identified in each of the genomes we are analyzing. The first step will be to run pileup which takes our multiple alignment files (.bam) and generates a VCF containing genotype likelihoods. The second step will be to make variant base calls. Some of the flags we utilize will be program specific while others will be our selections based upon the type of mutation predictions we have. The output file type will be a VCF file which stands for the variant call format.

[Link to bcftools](http://samtools.github.io/bcftools/bcftools.html)

Similar to before, we are going to construct a script containing 2 `for loops`.

~~~
#!/bin/bash
# what I do:

# Loop 1: mpileup on each sample
# Example of command for a single sample
# bcftools mpileup -Q 30 -C 50 -P Illumina \
# -Ov -B -f ref_genome/bdgp6.28.fa \
# -o vcfs/A44.pileup.vcf \
# clean_bams/A44.markdup.bam &

# -Q of 50 (minimum base quality for a base to be considered)
# -C of 50 (recommend value for BWA)

# Loop 2: sample call for each sample
# Example of a command for a single sample
# bcftools call -c -Ov -v -V 'indels' \
# -o vcfs/A44.calls.vcf vcfs/A44.pileup.vcf &

# -c is the original samtools/bcftools calling method
# -Ov is our output is an uncompressed vcf
# -v ouput variant sites only
# -V don't call indels = insertions or deletions

# Run me in clean_bams where the *.markdup.bam files are
~~~
{: .bash}


> ## Exercise
>
> Now that you have an example of the 2 `for loops` provided above, turn these into `for loops` that will run all of our samples. Name this new script
> `bcftools_mpileup.sh`. Be sure to include the output directory at the top of the script as follows:
> outdir=`../vcfs`  (The vcfs directory will be in results and where the output for the 2 `for loops` will be placed)
> mkdir -p $outdir
>
> When you are confident in your script, make it executable and run. this script will take a few hours as it is identifying SNPs in each genome.
>
>> ## Solution
>>
>> ~~~
>> cd ~/data/FlyCURE/scripts
>> nano bcftools_mpileup.sh
>> ~~~
>> {: .bash}
>>
>> ~~~
>> # Run me in clean_bams where the *.markdup.bam files are
>>
>> outdir='../vcfs'
>> mkdir -p $outdir
>>
>> # loop #1 bcftools mpileup
>> for i in *.markdup.bam; do
>>  prefix=$(basename $i .markdup.bam)
>>  bcftools mpileup -Q 30 -C 50 -P Illumina \
>>  -Ov -B -f ../../ref_genome/bdgp6.28.fa \
>>  -o $outdir/${prefix}.pileup.vcf \
>>  ${i} &
>> done
>> wait
>>
>> # loop #2 bcftools call
>> for i in *.markdup.bam; do
>>  prefix=$(basename $i .markdup.bam)
>>  bcftools call -c -Ov -v -V 'indels' \
>>  -o $outdir/${prefix}.calls.vcf $outdir/${prefix}.pileup.vcf &
>> done
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ chmod +x bcftools_mpileup.sh
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ cd ~/data/FlyCURE/results/clean_bams
>> $ ../../scripts/bcftools_mpileup.sh &
>> ~~~
>> {: .bash}
>>
> {: .solution}
{: .challenge}

This will likely take a few hours (4+) to run. I utilize `top` in one terminal to monitor that bcftools is running. I also utilize a second terminal to occasionally use `ls -lh` to see that files are being created and that the file sizes are changing. Once your output looks like the image below, you may proceed to reduce the file sizes, by zipping your `*.pileup.vcf`.

~~~
$ ~/data/FlyCURE/results/vcfs
$ ls -lh
~~~
{: .bash}

When the script is finished here is what you should see. In `top`, bcftools will no longer be present.

~~~
total 166G
-rw-r--r-- 1 gea_user gea_user   15 Dec  1 00:44 2R.bed
-rw-r--r-- 1 gea_user gea_user 110M Dec  1 22:59 A44.calls.vcf
-rw-r--r-- 1 gea_user gea_user  16G Dec  1 21:45 A44.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 119M Dec  1 22:59 B-2-13_S1.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:17 B-2-13_S1.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 116M Dec  1 22:58 B-2-16_S2.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:12 B-2-16_S2.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 117M Dec  1 22:59 Control.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 21:53 Control.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 132M Dec  1 23:00 cos2.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:19 cos2.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 128M Dec  1 22:59 H22.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:31 H22.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 125M Dec  1 22:59 L31.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:23 L31.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 125M Dec  1 22:59 L-3-2_S3.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:25 L-3-2_S3.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 121M Dec  1 22:59 N-1-1_S4.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:19 N-1-1_S4.pileup.vcf
-rw-r--r-- 1 gea_user gea_user 118M Dec  1 22:59 N-1-4_S5.calls.vcf
-rw-r--r-- 1 gea_user gea_user  17G Dec  1 22:23 N-1-4_S5.pileup.vcf
~~~
{: .output}

Reduce file size by zipping the `*.pileup.vcf`. Again, this will take some time. 
~~~
$ cd ~/data/FlyCURE/results/vcfs
$ gzip *.pileup.vcf
~~~
{: .bash}

**RUN CELL TWO OF YOUR PERSISTANCE NOTEBOOK and let it complete before closing.**
