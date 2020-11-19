---
layout: page
title: "Instructor Notes"
---

# Instructor Notes for Wrangling Genomics

## Issues with Macs vs Windows

Users are required to open *multiple* html files locally on their own browser and OS - will vary between users! Probably ctrl-click multiple selection within a file browser and right click will work on most systems?


## SAMtools or IGV?
Some instructors chose to use SAMtools tview for visualization of variant calling results while other prefer than IGV. SAMtools is the default because installation of IGV can take up additional instruction time, and SAMtools tview is sufficient to visualize results. However, episode 04-variant_calling includes instructions for installation and using IGV.

## Commands with Lengthy Run Times

#### Raw Data Downloads
The fastq files take about 15 minutes to download. This would be a good time to discuss the overall workflow of this lesson as illustrated by the graphic integrated on the page. It is recommended to start this lesson with the commmands to make and move to the /data/untrimmed-fastq directory and begin the download, and while files download, cover the "Bioinformatics Workflows" and "Starting with Data" texts. Beware that the last fastq file in the list takes the longest to download (~6-8 mins).

#### Running FastQC
The FastQC analysis on all raw reads takes about 10 minutes to run. It is a good idea to have learners start this command and cover the FastQC background material and images while FastQC runs.

#### Trimmomatic
The trimmomatic for loop will take about 10 minutes to run. Perhaps this would be a good time for a coffee break or a discussion about trimming.

#### bcftools mpileup
The bcftools mpileup command will take about 5 minutes to run. It is:

~~~
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
~~~

## Commands that must be modified
There are several commands that are example commands that will not run correctly if copy and pasted directly to the terminal. These commands serve as example commands and will need to be modified to fit each user. There is text around the commands outlining how they need to be changed, but it's helpful to be aware of them ahead of time as an instructor so you can set them up properly.

#### scp Command to Download FastQC to local machines
In the FastQC section, learners will download FastQC output files in order to open '.html `.html` summary files on their local machines in a web browser. The scp command currently contains a public DNS (for example, `ec2-34-238-162-94.compute-1.amazonaws.com`), but this will need to be replaced with the public DNS of the machine used by each learner. The Public DNS for each learner will be the same one they use to log in. The password will be provided to the Instructor when they receive instance information and will be the same for all learners.  

Command as is: 
~~~
scp dcuser@ec2-34-238-162-94.compute-1.amazonaws.com:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~

Command for learners to use:
~~~
scp dcuser@<Public_DNS>:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~

#### The unzip for loop
The for loop to unzip FastQC output will not work as directly copied pasted as:

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~

Because the `>` symbol will cause a syntax error when copied. This command will work correctly when typed at the command line! Learners may be surprised that a for loop takes multiple lines on the terminal.

#### unzip in Working with FastQC Output
The command `unzip *.zip` in the Working with FastQC Output section will run successfully for the first file, but fail for subsequent files. This error introduces the need for a for loop.

#### Example Trimmomatic Command
The first trimmomatic serves as an explanation for trimmomatic parameters and is not meant to be run. The command is:

~~~
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~

The correct syntax is outlined in the next section, Running Trimmomatic.

#### Actual Trimmomatic Command
The actual trimmomatic command is complicated for loop. It will need to be typed out by learners because the `>` symbols will raise an error if copy and pasted.

For reference, this command is:

~~~
$ for infile in *_1.fastq.gz
> do
>   base=$(basename ${infile} _1.fastq.gz)
>   trimmomatic PE ${infile} ${base}_2.fastq.gz \
>                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
>                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
>                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
> done
~~~

#### bwa mem Example Command
The first bwa mem command is an example and is not meant to be run. It is:

~~~
# bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sam
~~~

The correct command follows:

~~~
$ bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam
~~~

#### The Automation Episode
The code blocks at the beginning of the automation episode (05-automation.md) are examples of for loops and scripts and are not meant to be run by learners. The first code chunks that should be run are under Analyzing Quality with FastQC. 

Also, after the first code chunk of code meant to be run, there is a line that reads only `read_qc.sh` and will yield a message saying that this command wasn't found. After the creation of the script, this command will run the script that will be written.

