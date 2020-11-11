---
title: "Fly-CURE Assessing Read Quality"
teaching: 0
exercises: 0
questions:
- "How can I describe the quality of my data?"
objectives:
- "Explain how a FASTQ file encodes per-base quality scores."
- "Interpret a FastQC plot summarizing per-base quality across all reads."
- "Use `for` loops to automate operations on multiple files."
- "Run FastQC on fly NGS FASTQ files"
keypoints:
- "`for` loops let you perform the same set of operations on multiple files with a single command."
- "FastQC enables us to validate the continued use of the sequencing data"
---

# Bioinformatic workflows

When working with high-throughput sequencing data, the raw reads you get off of the sequencer will need to pass
through a number of  different tools in order to generate your final desired output. The execution of this set of
tools in a specified order is commonly referred to as a *workflow* or a *pipeline*.

An example of the workflow we will be using for our variant calling analysis is provided below with a brief
description of each step.

![workflow](../img/variant_calling_workflow.png)


1. Quality control - Assessing quality using FastQC
2. Quality control - Trimming and/or filtering reads (if necessary)
3. Align reads to reference genome
4. Perform post-alignment clean-up
5. Variant calling

These workflows in bioinformatics adopt a plug-and-play approach in that the output of one tool can be easily
used as input to another tool without any extensive configuration. Having standards for data formats is what
makes this feasible. Standards ensure that data is stored in a way that is generally accepted and agreed upon
within the community. The tools that are used to analyze data at different stages of the workflow are therefore
built under the assumption that the data will be provided in a specific format.  

# Starting with Data

Often times, the first step in a bioinformatic workflow is getting the data you want to work with onto a computer where you can work with it. If you have outsourced sequencing of your data, the sequencing center will usually provide you with a link that you can use to download your data. Today we will be working with private data shared with you through our Jupyter app as this is still unpublished data.

A raw_illumina directory was received directly from the Illumina sequencer. There are 4 quadrants on the flow cell and we want to merge these into one `.fastq` file to work with. This file merging has been done in advance for you by utilizing the cat and redirection command.

Launch your app as you always have with the dataset set to: `/iplant/home/your_username/data`. (Replace the `your_username` with your actual username). Once the app is launched, open and run the first cell of the Persistance Notebook and enter your password. Once that is completed, open the terminal and enter the following commands. Replace `your_username` with your actual username.

~~~
$ irysnc -r i:/iplant/home/kbieser/FlyCURE \
i:/iplant/home/your_username/data/FlyCURE
~~~
{: .bash}

Expect this to take a little while. To see that this was successful, if you return to the Discovery Environment home page and select `Data`, you should see a FlyCURE folder in your data folder.

Return to your Persistance Notebook and run cell 2. Once that cell completes running, `Complete and Save Outputs` for the app. Now you may re-launch the app using the same dataset as always (`/iplant/home/your_username/data`). (Continue to launch this way and run your Persistance Notebook each time). This will take a little while as you now have 10 new genomic sequences in your directory. Once JupyterLab opens you should now see a FlyCURE directory containing two subdirectories: `fastq_joined` (the paired end fastq sequencing files for each of the 10 fly stocks) and `adapters` (a file you will need for running FastQC).

The data in the `fastq_joined` directory comes in a compressed format, which is why there is a `.gz` at the end of the file names. This makes it faster to transfer, and allows it to take up less space on our computer. Let's unzip one of the files so that we can look at the fastq format.

If you recall, we can use the gunzip command to unzip '.gz' files. Since these files are much larger than the data sets we have worked with previously, it may take some time to unzip.  

~~~
$ cd ~/data/FlyCURE/fastq_joined
$ gunzip B-2-13_S1_R1.fastq.gz
~~~
{: .bash}

# Quality Control

We will now assess the quality of the sequence reads contained in our fastq files.

![workflow_qc](../img/var_calling_workflow_qc.png)
## Details on the FASTQ format

Although it looks complicated (and it is), we can understand FASTQ format with a little decoding. [Click here for a link to additional information about the FASTQ format.](https://en.wikipedia.org/wiki/FASTQ_format) Some rules about the format
include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

We can view the first complete read in one of the files our dataset by using `head` to look at
the first four lines.

~~~
$ head -n 4 B-2-13_S1_R1.fastq
~~~
{: .bash}

~~~
@NB501319:154:H332GAFX2:1:11101:23672:1019 1:N:0:CGATGT
CATTTTTTATATATTATATATTTTCTTGAAGCTGGTGCATATTCTATTTAGTCATCTTTTTTTATAGTGTCCTCAC
+
AAAAAEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEAE
~~~
{: .output}

Line 4 shows the quality for each nucleotide in the read. Quality is interpreted as the
probability of an incorrect base call (e.g. 1 in 10) or, equivalently, the base call
accuracy (e.g. 90%). To make it possible to line up each individual nucleotide with its quality
score, the numerical score is converted into a code where each individual character
represents the numerical quality score for an individual nucleotide. For example, in the line
above, the quality score line is:

~~~
AAAAAEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEAE
~~~
{: .output}

The numerical value assigned to each of these characters depends on the
sequencing platform that generated the reads. The sequencing machine used to generate our data
uses the standard Sanger quality PHRED score encoding, using Illumina version 1.8 onwards.
Each character is assigned a quality score between 0 and 41 as shown in
the chart below.

~~~
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
                   |         |         |         |         |
Quality score:    01........11........21........31........41                                
~~~
{: .output}

Each quality score represents the probability that the corresponding nucleotide call is
incorrect. This quality score is logarithmically based, so a quality score of 10 reflects a
base call accuracy of 90%, but a quality score of 20 reflects a base call accuracy of 99%.
These probability values are the results from the base calling algorithm and depend on how
much signal was captured for the base incorporation.

Looking back at our read:

~~~
@NB501319:154:H332GAFX2:1:11101:23672:1019 1:N:0:CGATGT
CATTTTTTATATATTATATATTTTCTTGAAGCTGGTGCATATTCTATTTAGTCATCTTTTTTTATAGTGTCCTCAC
+
AAAAAEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEAE
~~~
{: .output}

we can now see that there is a small range of quality scores.

> ## Exercise
>
> What is the last read in the `B-2-13_S1_R1.fastq` file? How confident
> are you in this read?
>
>> ## Solution
>> ~~~
>> $ tail -n 4 B-2-13_S1_R1.fastq
>> ~~~
>> {: .bash}
>>
>> ~~~
>> @NB501319:154:H332GAFX2:4:21612:17940:20399 1:N:0:GGATGT
>> ACTATGNATAAGCAAGAGAAAGCAAGCCCCAATGCAAATGGTTATAGATATCAGGAGCGCTCGGGGCGGGGGGGGC
>> +
>> AAAAAE#AEEEEEEEEEEEEEEE66EEE/E66/E//E6E6E//AAAEE////////////E////A/////</E/A
>> ~~~
>> {: .output}
>>
>> This read has the same consistent quality at the beginning as compared to the first read that we looked at, but has a wider range of quality scores, with scores becoming lower at the middle to end of the read. We will look at variations in position-based quality
>> in just a moment.
>>
> {: .solution}
{: .challenge}

At this point, lets validate that all the relevant tools are installed. If you are using JupyterLab then these _should_ be preinstalled.

~~~
$ fastqc -h
            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of
    which will help to identify a different potential type of problem in your
    data.

    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.

    The options for the program as as follows:

    -h --help       Print this help file and exit

    -v --version    Print the version of the program and exit

    -o --outdir     Create all output files in the specified output directory.                                                                                    
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.

    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.

    --nano          Files come from naopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.

    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.

    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created.  By default
                    this option will be set if fastqc is run in non-interactive
                    mode.

    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.

    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.

    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!

    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq

    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine

    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.

    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.

   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.

   -q --quiet       Supress all progress messages on stdout and only report errors.

   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.

BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
~~~
{: .bash}

if fastqc is not installed then you would expect to see an error like

~~~
$ fastqc -h
The program 'fastqc' is currently not installed. You can install it by typing:
sudo apt-get install fastqc
~~~
{: .bash}

If this happens check with your instructor before trying to install it.

## Assessing Quality using FastQC
In real life, you won't be assessing the quality of your reads by visually inspecting your
FASTQ files. Rather, you'll be using a software program to assess read quality and
filter out poor quality reads. We'll first use a program called FastQC to visualize the quality of our reads. [Click this link to take you to FastQC.](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Later in our workflow, we'll use another program to filter out poor quality reads.

FastQC has a number of features which can give you a quick impression of any problems your
data may have, so you can take these issues into consideration before moving forward with your
analyses. Rather than looking at quality scores for each individual read, FastQC looks at
quality collectively across all reads within a sample. The image below shows one FastQC-generated plot that indicates
a very high quality sample:

![good_quality](../img/good_quality1.8.png)

The x-axis displays the base position in the read, and the y-axis shows quality scores. In this
example, the sample contains reads that are 40 bp long. This is much shorter than the reads we
are working with in our workflow. For each position, there is a box-and-whisker plot showing
the distribution of quality scores for all reads at that position. The horizontal red line
indicates the median quality score and the yellow box shows the 1st to
3rd quartile range. This means that 50% of reads have a quality score that falls within the
range of the yellow box at that position. The whiskers show the absolute range, which covers
the lowest (0th quartile) to highest (4th quartile) values.

For each position in this sample, the quality values do not drop much lower than 32. This
is a high quality score. The plot background is also color-coded to identify good (green),
acceptable (yellow), and bad (red) quality scores.

Now let's take a look at a quality plot on the other end of the spectrum.

![bad_quality](../img/bad_quality1.8.png)

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the "bad" range, particularly on the tail end of the reads. The FastQC tool produces several other diagnostic plots to assess sample quality, in addition to the one plotted above.

## Running FastQC individual method (script method below - your choice)  

We will now assess the quality of the reads that have been shared with you. First, make sure you're still in the `fastq_joined` directory

~~~
$ cd ~/data/FlyCURE/fastq_joined
~~~
{: .bash}

> ## Exercise
>
>  How big are the files?
> (Hint: Look at the options for the `ls` command to see how to show
> file sizes.)
>
>> ## Solution
>>  
>> ~~~
>> $ ls -lh
>> ~~~
>> {: .bash}
>>
>> ~~~
>> -rw-r--r-- 1 gea_user gea_user 6.1G Oct  5 21:10 B-2-13_S1_R1.fastq
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:10 B-2-13_S1_R2.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.3G Oct  5 21:10 B-2-16_S2_R1.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.3G Oct  5 21:11 B-2-16_S2_R2.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:11 L-3-2_S3_R1.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:11 L-3-2_S3_R2.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:11 N-1-1_S4_R1.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:12 N-1-1_S4_R2.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.4G Oct  5 21:12 N-1-4_S5_R1.fastq.gz
>> -rw-r--r-- 1 gea_user gea_user 1.5G Oct  5 21:12 N-1-4_S5_R2.fastq.gz
>> ~~~
>> {: .output}
>>
>> There are ten FASTQ files ranging from 1.4G (1.4 gigabytes) to 6.1G.
>>
> {: .solution}
{: .challenge}

FastQC can accept multiple file names as input, and on both zipped and unzipped files, so we can use the *.fastq* wildcard to run FastQC on all of the FASTQ files in this directory.

~~~
$ fastqc *.fastq*
~~~
{: .bash}

You will see an automatically updating output message telling you the
progress of the analysis. It will start like this:

~~~
Started analysis of B-2-13_S1_R1.fastq
Approx 5% complete for B-2-13_S1_R1.fastq
Approx 10% complete for B-2-13_S1_R1.fastq
~~~
{: .output}

In total, this will take some time for FastQC to run on all
ten of our FASTQ files. When the analysis completes, your prompt
will return. So your screen will look something like this:

~~~
Approx 70% complete for N-1-4_S5_R2.fastq.gz
Approx 75% complete for N-1-4_S5_R2.fastq.gz
Approx 80% complete for N-1-4_S5_R2.fastq.gz
Approx 85% complete for N-1-4_S5_R2.fastq.gz
Approx 90% complete for N-1-4_S5_R2.fastq.gz
Approx 95% complete for N-1-4_S5_R2.fastq.gz
Analysis complete for N-1-4_S5_R2.fastq.gz
kbieser@aed3605cf.cyverse.run$
~~~
{: .output}

The FastQC program has created several new files within our
`~/fastq_joined` directory.

~~~
$ ls
~~~
{: .bash}

~~~
B-2-13_S1_R1.fastq        B-2-13_S1_R2_fastqc.zip   B-2-16_S2_R1.fastq.gz     L-3-2_S3_R1_fastqc.html  L-3-2_S3_R2_fastqc.zip   N-1-1_S4_R1.fastq.gz     N-1-4_S5_R1_fastqc.html  N-1-4_S5_R2_fastqc.zip
B-2-13_S1_R1_fastqc.html  B-2-13_S1_R2.fastq.gz     B-2-16_S2_R2_fastqc.html  L-3-2_S3_R1_fastqc.zip   L-3-2_S3_R2.fastq.gz     N-1-1_S4_R2_fastqc.html  N-1-4_S5_R1_fastqc.zip   N-1-4_S5_R2.fastq.gz
B-2-13_S1_R1_fastqc.zip   B-2-16_S2_R1_fastqc.html  B-2-16_S2_R2_fastqc.zip   L-3-2_S3_R1.fastq.gz     N-1-1_S4_R1_fastqc.html  N-1-1_S4_R2_fastqc.zip   N-1-4_S5_R1.fastq.gz
B-2-13_S1_R2_fastqc.html  B-2-16_S2_R1_fastqc.zip   B-2-16_S2_R2.fastq.gz     L-3-2_S3_R2_fastqc.html  N-1-1_S4_R1_fastqc.zip   N-1-1_S4_R2.fastq.gz     N-1-4_S5_R2_fastqc.html
~~~
{: .output}

For each input FASTQ file, FastQC has created a `.zip` file and a
`.html` file. The `.zip` file extension indicates that this is
actually a compressed set of multiple output files. We'll be working
with these output files soon. The `.html` file is a stable webpage
displaying the summary report for each of our samples.

We want to keep our data files and our results files separate, so we
will move these
output files into a new directory within our `results/` directory.

~~~
$ mkdir -p ~/data/Fly_CURE/results/fastqc_untrimmed_reads
$ mv *.zip ~/data/Fly_CURE/results/fastqc_untrimmed_reads/
$ mv *.html ~/data/Fly_CURE/results/fastqc_untrimmed_reads/
~~~
{: .bash}

Now we can navigate into this results directory and do some closer
inspection of our output files.

~~~
$ cd ~/data/Fly_CURE/results/fastqc_untrimmed_reads/
~~~
{: .bash}

## Running FastQC script method

If you want to run FastQC on all of your files at once and move the outputs to the correct directory, you can write a script to do so.

Let's start with making a scripts directory to save all the scripts we will be making

~~~
$ cd ~/data/FlyCURE
$ mkdir scripts
~~~
{: .bash}

Move into the scripts directory

~~~
$ cd scripts
~~~
{: .bash}

~~~
$ nano fastqc.sh
~~~
{: .bash}

In nano type in each command that you want to complete. Let's also make some notes for ourselves so that we know what to use this script for and where to run it from.

~~~
#!/bin/bash
# Use this script to run fastqc on my raw untrimmed reads
# Run me inside of the fastq_trimmed reads directory

mkdir -p ~/data/FlyCURE/results/fastqc_untrimmed_reads

fastqc -t 10 -o ../results/fastqc_untrimmed_reads *.fastq*

~~~
{: .bash}

In the FastQC command, the `-t 10`tells the server how many CPU's to utilize while running the command. The `-o ../results/fastqc_untrimmed_reads`, directs where the outputs from fastqc should be saved. In the past we have used the `mv` command after we ran fastqc, but this method allows the data to be placed in the correct directory upon creation rather than having to move them after the fact. It can do this because we prompted it to make the new directory first ensuring it exists.

Save the script.

To run the script you must make it executable.

~~~
$ chmod +x fastqc.sh
~~~
{: .bash}

Now run the script! Navigate to the correct directory where the script should be run and then call the script.

~~~
$ cd ~/data/FlyCURE/fastq_joined
$ ../scripts/fastqc.sh &
~~~
{: .bash}

The `&` will let you follow the progress. This will take some time to complete (mine took about 25 minutes) as there are 10 FASTQ files representing 1 ark control fly and 9 different Ethylmethanesulfonate (EMS) mutants. Remember these are the sequencing reads for the whole genome of _Drosophila melanogaster_ which is approximately 139.5 million base pairs.

Once completed, navigate into the results directory and do some closer
inspection of our output files.

~~~
$ cd ~/data/FlyCURE/results/fastqc_untrimmed_reads
~~~
{: .bash}

## Viewing the FastQC results

If we were working on our local computers, we'd be able to look at
each of these HTML files by opening them in a web browser.

However, these files are currently sitting on our remote CyVerse
instance, where our local computer can't see them.

The traditional way to look at these webpage summary reports would be
to transfer them to our local computers (i.e. your laptop). Since we are utilizing the JupyterLab app, we can just navigate to the directory and view them directly in JupyterLab.

In the left hand side of your console, navigate to the .html files.

~~~
data/FlyCURE/results/fastqc_untrimmed_reads
~~~

~~~
A44_R1_fastqc.html        B-2-13_S1_R2_fastqc.zip   Control_R2_fastqc.html  H22_R1_fastqc.zip   L-3-2_S3_R1_fastqc.html  N-1-1_S4_R2_fastqc.zip
A44_R1_fastqc.zip         B-2-16_S2_R1_fastqc.html  Control_R2_fastqc.zip   H22_R2_fastqc.html  L-3-2_S3_R1_fastqc.zip   N-1-4_S5_R1_fastqc.html
A44_R2_fastqc.html        B-2-16_S2_R1_fastqc.zip   cos2_R1_fastqc.html     H22_R2_fastqc.zip   L-3-2_S3_R2_fastqc.html  N-1-4_S5_R1_fastqc.zip
A44_R2_fastqc.zip         B-2-16_S2_R2_fastqc.html  cos2_R1_fastqc.zip      L31_R1_fastqc.html  L-3-2_S3_R2_fastqc.zip   N-1-4_S5_R2_fastqc.html
B-2-13_S1_R1_fastqc.html  B-2-16_S2_R2_fastqc.zip   cos2_R2_fastqc.html     L31_R1_fastqc.zip   N-1-1_S4_R1_fastqc.html  N-1-4_S5_R2_fastqc.zip
B-2-13_S1_R1_fastqc.zip   Control_R1_fastqc.html    cos2_R2_fastqc.zip      L31_R2_fastqc.html  N-1-1_S4_R1_fastqc.zip
B-2-13_S1_R2_fastqc.html  Control_R1_fastqc.zip     H22_R1_fastqc.html      L31_R2_fastqc.zip   N-1-1_S4_R2_fastqc.html
~~~
{: .output}

Double-click on one of the .html files and a new tab will automatically open with the FastQC Report. Open and review each of the 10 fastqc.html files. View the next section for information about decoding the FastQC outputs.

> ## Exercise
>
> View each of the FastQC reports. Which sample(s) looks the best
> in terms of per base sequence quality? Which sample(s) look the
> worst? Write some notes about what you are seeing in your notebook.
>
>> ## Solution
>> Since we conducted paired-end sequencing, there are two reads for each fastq. A R1 (read 1) and a R2 (read 2). Each sequencing file has usable data, but the quality decreases >> toward the end of the reads. R2 has poorer quality than R1.
> {: .solution}
{: .challenge}

## Decoding the other FastQC outputs
We've now looked at quite a few "Per base sequence quality" FastQC graphs, but there are nine other graphs that we haven't talked about! Below we have provided a brief overview of interpretations for each of these plots. For more information, please see the FastQC documentation [by clicking this link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

+ [**Per tile sequence quality**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html): the machines that perform sequencing are divided into tiles. This plot displays patterns in base quality along these tiles. Consistently low scores are often found around the edges, but hot spots can also occur in the middle if an air bubble was introduced at some point during the run.
+ [**Per sequence quality scores**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html): a density plot of quality for all reads at all positions. This plot shows what quality scores are most common.
+ [**Per base sequence content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html): plots the proportion of each base position over all of the reads. Typically, we expect to see each base roughly 25% of the time at each position, but this often fails at the beginning or end of the read due to quality or adapter content.
+ [**Per sequence GC content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html): a density plot of average GC content in each of the reads.  
+ [**Per base N content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html): the percent of times that 'N' occurs at a position in all reads. If there is an increase at a particular position, this might indicate that something went wrong during sequencing.  
+ [**Sequence Length Distribution**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html): the distribution of sequence lengths of all reads in the file. If the data is raw, there is often on sharp peak, however if the reads have been trimmed, there may be a distribution of shorter lengths.
+ [**Sequence Duplication Levels**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html): A distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.  
+ [**Overrepresented sequences**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html): A list of sequences that occur more frequently than would be expected by chance.
+ [**Adapter Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html): a graph indicating where adapater sequences occur in the reads.
+ [**K-mer Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html): a graph showing any sequences which may show a positional bias within the reads.
