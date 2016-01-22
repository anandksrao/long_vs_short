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
The aim of this exercise is to demonstrate the advantages of long reads in the assembly of highly repetitive genome sequences. We will try to assemble a small part of the chimp genome known to be highly repetitive.

### Generate test reads
For the purposes of this demonstration we will use simulated data. If you would like to skip this step you can find a dataset for long and short reads here.

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

Perform *de-novo* assembly with long reads using Canu.

```
canu -p longReadsCov30 -d longReadsCov30 genomesize=184664 -s specfile.dat -nanopore-raw longReadsCov30.fasta
```

### Assessment of assemblies
Assess both of your assemblies with QUAST:

```
quast.py consensusFile/contigFile -R refChimp1.fna
```

Open quast_results/results/report.pdf to see the results. What do you see?
Spoiler: You should see a large difference in '% genome covered'.

### Optional: How does alignment of short reads compare to *de-novo* assembly of short reads in terms of '% genome covered'?

To align the short reads we need a different alignment tool. Here, we will use Bowtie 2.

```
bowtie2-build reference_chimp1.fna reference_chimp1.fna
bowtie2 -f -x reference_chimp1.fna -U shortReadsCov30.fasta -S shortReadsCov30_aligned.sam
samtools view -bS shortReadsCov30_aligned.sam > shortReadsCov30_aligned.bam
samtools sort shortReadsCov30_aligned.bam shortReadsCov30_aligned.sorted
samtools index shortReadsCov30_aligned.sorted.bam
```

To see how much of your genome was mapped, run:

```
samtools flagstat shortReadsCov30_aligned.sorted.bam
```

How do you explain the difference in '% genome covered' between the *de-novo* and the alignment assembly?

