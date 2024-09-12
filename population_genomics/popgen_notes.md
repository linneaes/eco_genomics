# Coding and data notes for the population genomics module

## Author: Linnea Ericsson Slater

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCFs)

We learned how to set up our notes file and use markdown shortcuts. We then used the command line by opening VACC Shell Access and learned how to switch our directories over to the example data, then how to open it with zcat and use \| head to look at the first 10 lines. We then learned how to change the amount of data lines in the head by adding -n and a number (example: zcat file.gz \| head -n 4).

cd = change directory

zcat is used to view a zip file

\|s -\| = long list -\> can be abbreviated as ll

head = view first 10 lines of script

\| "pipe" = send output of one function to another

If you are in the correct directory w/ a specific file or directory, you can use tab to auto-complete the name

Our example data was a Q-file, so we saw the Q-scores in I's and G's according to what the score was, along with the sequence itself, and the barcode identifying the individual sample.

### 09-12-2024 - Viewing VCF files and talking about filtering

We learned how to open samtools on the command line and use it to look at a bam file to see the summary of reads and the sequencing depth from the data.

cd .. gets you out of a directory one time back

q to quit
