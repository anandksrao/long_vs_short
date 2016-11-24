# Comparing assembly of repetitive sequences using long and short reads.

**Software requirements**
- [ReadSim](http://sourceforge.net/p/readsim/wiki/manual/) (read simulator)
- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) (short read *de-novo* assembler)
- [Canu](https://github.com/marbl/canu/releases) (long read *de-novo* assembler)
- [QUAST](http://bioinf.spbau.ru/quast) (assessment of assemblies)

**Optional:**
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (short read aligner)
- [LAST](http://last.cbrc.jp/) (long read aligner)
- [samtools](http://www.htslib.org/download/) (for manipulating alignments)


### Introduction
The aim of this exercise is to demonstrate the advantages of long reads in the assembly of difficult genome sequences. We will try to assemble a small part of the chimp genome known to be highly repetitive. For a tutorial on de-novo assembly and read alignment see [here](https://github.com/demharters/assemblyTutorial).

### Generate test reads
For the purposes of this demonstration we will use simulated data. Download the reference sequence from [here](https://figshare.com/s/be47dc169f8759545b5a).
*If you would like to skip this step you can find simulated datasets for long and short reads under the same link.*

Generate short reads with the following command:

```
./readsim.py sim fa --ref refChimp1.fna --pre shortReadsCov30 --rev_strd on
  --tech nanopore --read_mu 30 --read_dist normal --cov_mu 30
  --err_sub_mu 0.001 --err_in_mu 0.001 --err_del_mu 0.001
```
*Note these are double dashes "--"*

This simulation will generate a set of short fasta reads (30 bases on average) with a 30x coverage using our reference sequence as a template. We set the substitution, insertion and deletion error rates to 0.1% to replicate the typical characteristics of short reads (even though the technology is set to "nanopore").

Generate long reads with the following command:

```
/readsim.py sim fa --ref refChimp1.fna --pre longReadsCov30
  --rev_strd on --tech nanopore --read_mu 15000 --read_dist normal --cov_mu 30
```

This simulation will generate a set of long fasta reads (15kb on average) also with a 30x coverage and using our reference sequence as a template. The substitution, insertion and deletion error rates are set to 3% (i.e. we are simulating corrected reads). The error rates selected here are high and are constantly improving as Nanopore technology matures. For more information on correction see (Loman 2015, Nature Methods)[http://www.nature.com/nmeth/journal/v12/n8/full/nmeth.3444.html].


### *De-novo* assembly

Perform *de-novo* assembly with short reads using Velvet.

``` 
velveth shortReadsCov30_assembly 21 shortReadsCov30.fasta
velvetg shorReadsCov30_assembly
```
This will create a folder called “shortReadsCov30_assembly” that contains your short-read assembly.

Velveth takes a number of sequence files as input, generates a hashtable and spits out two files sequences and roadmaps, which are required by velvetg.

In the velveth command '21' is the kmer length used for the hash table. If you would like to know more read section 5.2 in the [manual](http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf). If you are interested try playing around with this value and see how it affects the assembly.

For more information on how velvet works, see [here](http://microbialinformaticsj.biomedcentral.com/articles/10.1186/2042-5783-3-2).

Perform *de-novo* assembly with long reads using Canu.

```
canu -p longReadsCov30 -d longReadsCov30_assembly genomesize=184664 -s specfile.dat -nanopore-raw longReadsCov30.fasta
```
*Canu requires a file "specfile.dat" in the working directory. This file is used to pass more options to canu. For our purposes this file can be empty.*

This will create a folder called “longReadsCov30_assembly” that contains your long-read assembly.

Canu is based on the Celera assembler. For terminology see [here](http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology).

### Assessment of assemblies
Assess both of your assemblies with QUAST:

```
quast.py assemblyFolder/contigFile -R refChimp1.fna
```
*The contig file will be named something like "longReadsCov30.contig.fa" or "contigs.fa".*

Open quast_results/results/report.pdf to see the results. What do you see?

*Spoiler: You should see a large difference in '% genome covered'.*


### Questions
- Q1. What is the main difference between the two assemblies?
Let's assess our reference sequence with RepeatMasker. This tool will allow us to identify any repetitive elements.
- Go to the [RepeatMasker webserver](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker) and upload the reference sequence.
*Note, the sequence has to be shorter than 100kb, so you will have to modify the file. If the run takes too long, try a shorter sequence.*
- Press 'submit sequence'.
As you can see our sequence is full of repetitive elements. This is bad news for short reads and it is most likely the reason for our poor *de-novo* assembly. Our long read assembly did just fine it seems. It demonstrates that assembly of repetitive genome regions requires long reads that span entire repeats.

- Q2. How do the results explain the difference in quality between the two assemblies?
- Q3. Can you identify the minimum read length required to assemble this chimp sequence by 80%? Assume a 30x coverage and 0.03 error rates. (Go back to the readsim.py step and simulate new reads).
- Q4a. How high can you go with the error rates in the simulated long reads before the assembly starts to fail? Assume 15kb read length and 30x coverage.
- Q4b. Once you managed to break the long read de-novo assembly by increasing the error rate, does tweaking one of the other parameters restore the assembly?
- Q5. How does alignment of short reads to a reference sequence compare to short read de-novo assembly in terms of '% genome covered'?

To assemble the short reads by alignment we need a different alignment tool. Here, we will use Bowtie 2.

```
bowtie2-build refChimp1.fna refChimp1.fna
bowtie2 -f -x refChimp1.fna -U shortReadsCov30.fasta -S shortReadsCov30_aligned.sam
```
<!--- samtools view -bS shortReadsCov30_aligned.sam > shortReadsCov30_aligned.bam --->
<!--- samtools sort shortReadsCov30_aligned.bam shortReadsCov30_aligned.sorted.bam --->
<!--- samtools index shortReadsCov30_aligned.sorted.bam --->
How much of the sequence has been covered (it should say in the terminal)?

<!--- Tutorial on [samtools](http://biobits.org/samtools_primer.html). --->

<!--- To see how much of your genome was mapped, run:

<!--- samtools flagstat shortReadsCov30_aligned.sorted.bam --->

How do you explain the difference in '% genome covered' between the *de-novo* and the alignment assembly?

