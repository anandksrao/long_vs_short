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
./readsim.py sim fa --ref refChimp1.fna --pre shortReadsCov30 \
--rev_strd on --tech nanopore --read_mu 30 --read_dist normal --cov_mu 30 \
--err_sub_mu 0.001 --err_in_mu 0.001 --err_del_mu 0.001
```

This simulation will generate a set of short fasta reads (30 bases on average) with a 30x coverage using our reference sequence as a template. We set the substitution, insertion and deletion error rates to 0.1% to replicate the typical characteristics of short reads (even though the technology is set to "nanopore").

Generate long reads with the following command:

```
/readsim.py sim fa --ref refChimp1.fna --pre longReadsCov30 \ 
--rev_strd on --tech nanopore --read_mu 15000 --read_dist normal --cov_mu 30
```

This simulation will generate a set of long fasta reads (15kb on average) also with a 30x coverage and using our reference sequence as a template. The substitution, insertion and deletion error rates are set to 3% (i.e. we are simulating corrected reads).


### *De-novo* assembly

Perform *de-novo* assembly with short reads using Velvet.

``` 
velveth chimp1 21 shortReadsCov30.fasta
velvetg chimp1
```
Velveth takes a number of sequence files as input, generates a hashtable and spits out two files sequences and roadmaps, which are required by velvetg.

In the velveth command '21' is the kmer length used for the hash table. If you would like to know more read section 5.2 in the [manual](http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf). If you are interested try playing around with this value and see how it affects the assembly.

For more information on how velvet works, see [here](http://microbialinformaticsj.biomedcentral.com/articles/10.1186/2042-5783-3-2).

Perform *de-novo* assembly with long reads using Canu.

```
canu -p longReadsCov30 -d longReadsCov30 genomesize=184664 -s specfile.dat -nanopore-raw longReadsCov30.fasta
```
Canu is based on the Celera assembler. For terminology see [here](http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology).

### Assessment of assemblies
Assess both of your assemblies with QUAST:

```
quast.py consensusFile/contigFile -R refChimp1.fna
```

Open quast_results/results/report.pdf to see the results. What do you see?

*Spoiler: You should see a large difference in '% genome covered'.*

### Why do the short reads fail?
Let's assess our genomic sequence with RepeatMasker. This tool will allow us to identify any repetitive elements.
Go to the [RepeatMasker webserver](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker) and upload your file. Press 'submit sequence'.
As you can see our sequence is full of repetitive elements. This is bad news for short reads and it is most likely the reason for our poor *de-novo* assembly. Our long read assembly did just fine it seems. It demonstrates that assembly of repetitive genome regions requires long reads that span the entire repeats.

### Questions
- Can you identify the minimum read length required to assemble this chimp sequence by 90%, 95%, 99%? Assume a 30x coverage and 0.001 error rates.
- How high can you go with the error rates in the simulated long reads before assembly starts to fail? Assume 15kb read length and 30x coverage. Once you managed to break the long read *de-novo* assembly, does increasing the coverage restore the assembly?

### Optional: How does alignment of short reads compare to *de-novo* assembly of short reads in terms of '% genome covered'?

To assemble the short reads by alignment we need a different alignment tool. Here, we will use Bowtie 2.

```
bowtie2-build reference_chimp1.fna reference_chimp1.fna
bowtie2 -f -x reference_chimp1.fna -U shortReadsCov30.fasta -S shortReadsCov30_aligned.sam
samtools view -bS shortReadsCov30_aligned.sam > shortReadsCov30_aligned.bam
samtools sort shortReadsCov30_aligned.bam -T shortReadsCov30_aligned.sorted -o shortReadsCov30_aligned.sorted.bam
samtools index shortReadsCov30_aligned.sorted.bam
```

Tutorial on [samtools](http://biobits.org/samtools_primer.html).

To see how much of your genome was mapped, run:

```
samtools flagstat shortReadsCov30_aligned.sorted.bam
```

How do you explain the difference in '% genome covered' between the *de-novo* and the alignment assembly?

