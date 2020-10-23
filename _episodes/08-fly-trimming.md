---
title: "Trimming and Filtering"
teaching: 0
exercises: 0
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using Trimmomatic."
- "Select and set multiple options for command-line bioinformatic tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---

# Cleaning Reads

In the previous lesson, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic Options

Trimmomatic has a variety of options to trim your reads. If we run the following command, we can see some of our options.

~~~
$ trimmomatic
~~~
{: .bash}

Which will give you the following output:
~~~
Usage:
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or:
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or:
       -version
~~~
{: .output}

This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads.
Next, we specify what flag we would like to run. For example, you can specify `threads` to indicate the number of
processors on your computer that you want Trimmomatic to use. In most cases using multiple threads (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important.
In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

| option    | meaning |
| ------- | ---------- |
|  \<inputFile1>  | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.|
| \<inputFile2> | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.|
|  \<outputFile1P> | Output file that contains surviving pairs from the `_1` file. |
|  \<outputFile1U> | Output file that contains orphaned reads from the `_1` file. |
|  \<outputFile2P> | Output file that contains surviving pairs from the `_2` file.|
|  \<outputFile2U> | Output file that contains orphaned reads from the `_2` file.|

The last thing trimmomatic expects to see is the trimming parameters:

| step   | meaning |
| ------- | ---------- |
| `ILLUMINACLIP` | Perform adapter removal. |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

However, a complete command for Trimmomatic will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

~~~
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `PE` | that it will be taking a paired end file as input |
| `-threads 4` | to use four computing threads to run (this will speed up our run) |
| `SRR_1056_1.fastq` | the first input file name |
| `SRR_1056_2.fastq` | the second input file name |
| `SRR_1056_1.trimmed.fastq` | the output file for surviving pairs from the `_1` file |
| `SRR_1056_1un.trimmed.fastq` | the output file for orphaned reads from the `_1` file |
| `SRR_1056_2.trimmed.fastq` | the output file for surviving pairs from the `_2` file |
| `SRR_1056_2un.trimmed.fastq` | the output file for orphaned reads from the `_2` file |
| `ILLUMINACLIP:SRR_adapters.fa`| to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
|`SLIDINGWINDOW:4:20` | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |



> ## Multi-line commands
> Some of the commands we ran in this lesson are long! When typing a long
> command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
{: .callout}



## Running Trimmomatic for a single sample

Now we will run Trimmomatic on our data. To begin, navigate to your `fastq_joined` data directory:

~~~
$ cd ~/data/FlyCURE/fastq_joined
~~~
{: .bash}

We are going to run Trimmomatic on one of our paired-end samples.
While using FastQC we are going to use a custom set of adapters to include the polyG that was present in some of our samples as an artifact. The rest of the adapaters in this file include NexteraPE, TruSeq2, and TruSeq3 adapter sequences that come standard with the installation of Trimmomatic. I have added this file to your `FlyCURE/adapters` directory and it is called `custom_PE.fa`.

We will use a leading and trailing threshold quality of 3. The entire read will be dropped if the sequence length is less than 36 bp (the MINLEN setting). This command will take about 10 minutes to run.

~~~
$ trimmomatic PE -threads 8 \
A44_R1.fastq.gz A44_R2.fastq.gz \
A44_R1.trim.fastq.gz A44_R1.orphan.fastq.gz \
A44_R2.trim.fastq.gz A44_R2.orphan.fastq.gz \
ILLUMINACLIP:../adapters/custom_PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
~~~
{: .bash}


~~~
TrimmomaticPE: Started with arguments:
TrimmomaticPE: Started with arguments:
 -threads 8 fastq_joined/A44_R1.fastq.gz fastq_joined/A44_R2.fastq.gz fastq_trimmed/A44_R1.trim.fastq.gz orphans/A44_R1.orphan.fastq.gz fastq_trimmed/A44_R2.trim.fastq.gz orphans/A44_R2.orphan.fastq.gz ILLUMINACLIP:../adapters/custom_PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'
Using Long Clipping Sequence: 'TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC'
Using Long Clipping Sequence: 'TTTTTTTTTTCAAGCAGAAGACGGCATACGA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
Using Long Clipping Sequence: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLUMINACLIP: Using 2 prefix pairs, 15 forward/reverse sequences, 1 forward only sequences, 1 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 19671500 Both Surviving: 18935218 (96.26%) Forward Only Surviving: 1 (0.00%) Reverse Only Surviving: 0 (0.00%) Dropped: 736281 (3.74%)
TrimmomaticPE: Completed successfully
~~~
{: .output}

You may have noticed that Trimmomatic automatically detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output files:

~~~
$ ls A44*
~~~
{: .bash}

~~~
A44_R1.fastq.gz  A44_R1.orphan.fastq.gz  A44_R1.trim.fastq.gz  A44_R2.fastq.gz  A44_R2.orphan.fastq.gz  A44_R2.trim.fastq.gz
~~~
{: .output}

The output files are also FASTQ files. The output .trim files should be smaller than our
input file, because we've removed reads. The output orphan files are the smallest because those are sequences that we removed because their mate was a terrible read. The absolute worst reads are now gone. It's partner ended up in those orphan files. We can confirm this:

~~~
$ ls A44* -lh
~~~
{: .bash}

~~~
-rw-r--r-- 1 gea_user gea_user 935M Oct 23 16:40 A44_R1.fastq.gz
-rw-r--r-- 1 gea_user gea_user  116 Oct 23 17:08 A44_R1.orphan.fastq.gz
-rw-r--r-- 1 gea_user gea_user 908M Oct 23 17:08 A44_R1.trim.fastq.gz
-rw-r--r-- 1 gea_user gea_user 956M Oct 23 16:40 A44_R2.fastq.gz
-rw-r--r-- 1 gea_user gea_user   20 Oct 23 17:08 A44_R2.orphan.fastq.gz
-rw-r--r-- 1 gea_user gea_user 908M Oct 23 17:08 A44_R2.trim.fastq.gz
~~~
{: .output}


We've just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly!

We unzipped one of our files before to work with it, let's compress it again before we run our for loop.

~~~
$ gzip B-2-13_S1_R1.fastq
~~~
{: .bash}

~~~
$ for infile in *_R1.fastq.gz
> do
> echo ${infile}
> done
~~~
{: .bash}

~~~
A44_R1.fastq.gz
B-2-13_S1_R1.fastq.gz
B-2-16_S2_R1.fastq.gz
Control_R1.fastq.gz
cos2_R1.fastq.gz
H22_R1.fastq.gz
L31_R1.fastq.gz
L-3-2_S3_R1.fastq.gz
N-1-1_S4_R1.fastq.gz
N-1-4_S5_R1.fastq.gz
~~~
{: .output}

~~~
$ for infile in *_R1.fastq.gz
> do
> echo ${infile}
> echo "this is a line that isn't a filename"
> done
~~~
{: .bash}

~~~
A44_R1.fastq.gz
this is a line that isn't a filename
B-2-13_S1_R1.fastq.gz
this is a line that isn't a filename
B-2-16_S2_R1.fastq.gz
this is a line that isn't a filename
Control_R1.fastq.gz
this is a line that isn't a filename
cos2_R1.fastq.gz
this is a line that isn't a filename
H22_R1.fastq.gz
this is a line that isn't a filename
L31_R1.fastq.gz
this is a line that isn't a filename
L-3-2_S3_R1.fastq.gz
this is a line that isn't a filename
N-1-1_S4_R1.fastq.gz
this is a line that isn't a filename
N-1-4_S5_R1.fastq.gz
this is a line that isn't a filename
~~~
{: .output}

Consider the original trimmomatic command, it had 2 inputs and 4 output files and our loop currently is only aware of 1 of our 2 input files. We need to think of a way to make it aware of all of the others. What do the file names have in common?

~~~
$ ls -lh
~~~
{: .bash}

It looks like there is a common suffix or base to the name. So we can manipulate the content of the variable that we already know called infile. We can manipulate the variable called infile and save those manipulated contents to a new variable.

Let's use basename to demonstrate that we can remove the suffix.

~~~
$ basename A44_R1.fastq.gz _R1.fastq.gz
~~~
{: .bash}

~~~
A44
~~~
{: .output}

How can we get the output of basename to get stored in a new variable? We can use a special syntax with $() to start a command, run something, and then get back to finishing what it started. Let's make a new variable and do a subshell of it to see how it works. The output of what is run in () gets stored in prefix.

~~~
$ prefix=$(basename A44_R1.fastq.gz _R1.fastq.gz)
$ echo ${prefix}
~~~
{: .bash}

~~~
A44
~~~
{: .output}

We can combine the basename command with our existing variable. Let's go to the script we are writing.

~~~
$ for infile in *_R1.fastq.gz
> do
>  echo ${infile}
>  base=$(basename ${infile} _R1.fastq.gz)
>  echo ${base}
> done
~~~
{: .bash}

~~~
A44_R1.fastq.gz
A44
B-2-13_S1_R1.fastq.gz
B-2-13_S1
B-2-16_S2_R1.fastq.gz
B-2-16_S2
Control_R1.fastq.gz
Control
cos2_R1.fastq.gz
cos2
H22_R1.fastq.gz
H22
L31_R1.fastq.gz
L31
L-3-2_S3_R1.fastq.gz
L-3-2_S3
N-1-1_S4_R1.fastq.gz
N-1-1_S4
N-1-4_S5_R1.fastq.gz
N-1-4_S5
~~~
{: .output}

Remember your history command? We had a pretty complicated trimmomatic command that we used and we want to substitute the variable in there. Let's remember what that command was. Who remembers what command shows us things we've already done? Look for the trimmomatic command in their history.

~~~
$ history | grep trimmomatic
~~~
{: .bash}

Copy your trimmomatic command and open nano.

~~~
$ nano
~~~
{: .bash}

Nano is unfortunately bad with copying and pasting. You will need to scroll by using arrow keys and use a \ return to get the lines in a readable format.

~~~
trimmomatic PE -threads 8 \
  A44_R1.fastq.gz A44_R2.fastq.gz \
  A44_R1.trim.fastq.gz A44_R1.orphan.fastq.gz \
  A44_R2.trim.fastq.gz A44_R2.orphan.fastq.gz \
  ILLUMINACLIP:../adapters/custom_PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
~~~
{: .bash}

Once it is copied in you are now going to convert it to a for loop by replacing the explicit parts with variables.

> ## Exercise
> Using the steps we tested above, create a for loop to run trimmomatic on all of our fastq files. Since you have the trimmomatic command pasted into nano, comment out each line > and use it to help you build the for loop. Save the text file as trim.sh and move to your scripts directory. If you really want to challenge yourself, you can add in commands  > to make a new fastq_trimmed/ and move the .trim and .orphan files to that new directory all within the script. Comment out notes within the text to help yourself. Don't forget > you have to make the script executable in order to run > it. Run the for loop.  
>
>
>> ## Solution
>> I moved into scripts/ and used nano to open a trim.sh
>> ~~~
>> for infile in *_R1.fastq.gz
>> do
>> base=$(basename ${infile} _R1.fastq.gz)
>> trimmomatic PE -threads 8 \
>> ${infile} ${base}_R2.fastq.gz \
>> ${base}_R1.trim.fastq.gz ${base}_R1.orphan.fastq.gz \
>> ${base}_R2.trim.fastq.gz ${base}_R2.orphan.fastq.gz \
>> ILLUMINACLIP:../adapters/custom_PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
>> done
>>
>> mkdir -p ~/data/FlyCURE/results/fastq_trimmed
>> mv *.trim* ../results/fastq_trimmed
>> mv *.orphan* ../results/fastq_trimmed
>> ~~~
>> {: .bash}
>>
>> If you didn't create the script within the scripts directory, then:
>> ~~~
>> $ mv trim.sh ~/data/FlyCURE/scripts
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ chmod +x trim.sh
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ cd ~/data/FlyCURE/fastq_joined
>> $ ../scripts/trim.sh &
>> ~~~
>> {: .bash}
>> This took about 2.5 hours for me.
>> ~~~
>> Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
>> Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
>> Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
>> Using Long Clipping Sequence: 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
>> Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
>> Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
>> Using Long Clipping Sequence: 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'
>> Using Long Clipping Sequence: 'TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC'
>> Using Long Clipping Sequence: 'TTTTTTTTTTCAAGCAGAAGACGGCATACGA'
>> Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
>> Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
>> Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
>> Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
>> Using Long Clipping Sequence: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
>> ILLUMINACLIP: Using 2 prefix pairs, 15 forward/reverse sequences, 1 forward only sequences, 1 reverse only sequences
>> Quality encoding detected as phred33
>> Input Read Pairs: 19671500 Both Surviving: 18658198 (94.85%) Forward Only Surviving: 276410 (1.41%) Reverse Only Surviving: 39388 (0.20%) Dropped: 697504 (3.55%)
>> TrimmomaticPE: Completed successfully
>> ~~~
>> {: .output}
>>
> {: .solution}
{: .challenge}

~~~
$ cd ~/data/FlyCURE/results/fastq_trimmed
~~~
{: .bash}

~~~
A44_R1.orphan.fastq.gz        B-2-13_S1_R2.orphan.fastq.gz  Control_R1.orphan.fastq.gz  cos2_R2.orphan.fastq.gz  L31_R1.orphan.fastq.gz       L-3-2_S3_R2.orphan.fastq.gz  N-1-4_S5_R1.orphan.fastq.gz
A44_R1.trim.fastq.gz          B-2-13_S1_R2.trim.fastq.gz    Control_R1.trim.fastq.gz    cos2_R2.trim.fastq.gz    L31_R1.trim.fastq.gz         L-3-2_S3_R2.trim.fastq.gz    N-1-4_S5_R1.trim.fastq.gz
A44_R2.orphan.fastq.gz        B-2-16_S2_R1.orphan.fastq.gz  Control_R2.orphan.fastq.gz  H22_R1.orphan.fastq.gz   L31_R2.orphan.fastq.gz       N-1-1_S4_R1.orphan.fastq.gz  N-1-4_S5_R2.orphan.fastq.gz
A44_R2.trim.fastq.gz          B-2-16_S2_R1.trim.fastq.gz    Control_R2.trim.fastq.gz    H22_R1.trim.fastq.gz     L31_R2.trim.fastq.gz         N-1-1_S4_R1.trim.fastq.gz    N-1-4_S5_R2.trim.fastq.gz
B-2-13_S1_R1.orphan.fastq.gz  B-2-16_S2_R2.orphan.fastq.gz  cos2_R1.orphan.fastq.gz     H22_R2.orphan.fastq.gz   L-3-2_S3_R1.orphan.fastq.gz  N-1-1_S4_R2.orphan.fastq.gz
B-2-13_S1_R1.trim.fastq.gz    B-2-16_S2_R2.trim.fastq.gz    cos2_R1.trim.fastq.gz       H22_R2.trim.fastq.gz     L-3-2_S3_R1.trim.fastq.gz    N-1-1_S4_R2.trim.fastq.gz
~~~
{: .output}


> ## Exercise
>
> Now that our samples have gone through quality control, they should perform
> better on the quality tests run by FastQC. Go ahead and re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming. Modify your fastqc.sh script to run on the trim_fastq. Again, plan for this to take a few hours.
>
>> ## Solution
>>
>> First edit your fastqc.sh and save as fastqc_trimmed.sh in your scripts/.
>>
>> ~~~
>> fastqc *.trim.fastq*
>> mkdir -p ~/data/FlyCURE/results/fastqc_trimmed_reads
>> mv *.zip ~/data/FlyCURE/results/fastqc_trimmed_reads
>> mv *.html ~/data/FlyCURE/results/fastqc_trimmed_reads
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ chmod +x fastqc_trimmed.sh
>> ~~~
>> {: .bash}
>>
>> ~~~
>> $ cd ~/data/FlyCURE/results/fastq_trimmed
>> $ ../../scripts/fastqc_trimmed.sh &
>> ~~~
>> {: .bash}
>>
> {: .solution}
{: .challenge}
