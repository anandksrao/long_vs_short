# Comparing assembly of repetitive sequences using long and short reads.

Software requirements
- [ReadSim](http://sourceforge.net/p/readsim/wiki/manual/) (read simulator)
- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) (short read *de-novo* assembler)
- [Canu](https://github.com/marbl/canu/releases) (long read *de-novo* assembler)

  optional:
  - [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (short read aligner)
  - [LAST](http://last.cbrc.jp/) (long read aligner)
  - [samtools](http://www.htslib.org/download/) (for manipulating alignments)


### Introduction

### Generate test reads
For the purposes of this demonstration we will use simulated data. If you would like to skip this step you can find a dataset for long and short reads here.

Generate short reads with the following command:
```
./readsim.py sim fa --ref refChimp3.fna --pre chimp3_shortReadsCov30_2 \
--rev_strd on --tech nanopore --read_mu 30 --read_dist normal --cov_mu 30 \
--err_sub_mu 0.001 --err_in_mu 0.001 --err_del_mu 0.001
```

This simulation will generate a set of short reads (30 bases on average) with a 30x coverage using our reference sequence as a template. We set the substitution, insertion and deletion error rates to 0.1% to replicate the typical characteristics of short reads (even though the technology is set to "nanopore").

Generate long reads with the following command:

'''
/readsim.py sim fa --ref refChimp3.fna --pre chimp3_longReadsCov30 \ 
--rev_strd on --tech nanopore --read_mu 15000 --read_dist normal --cov_mu 30
'''

This simulation will generate a set of long reads (15kb on average) also with a 30x coverage and using our reference sequence as a template. The substitution, insertion and deletion error rates are set to 3% (i.e. we are simulating corrected reads).


### *De-novo* assembly


### Optional: Does alignment to a reference change things?

To align the short reads we need a different alignment tool. Here, we will use Bowtie 2.

```
bowtie2-build reference_chimp1.fna reference_chimp1.fna
bowtie2 -f -x reference_chimp1.fna -U shortReads_chimp1.fasta -S shortReads_chimp1_aligned.sam
samtools view -bS shortReads_chimp1_aligned.sam > shortReads_chimp1_aligned.bam
samtools sort shortReads_chimp1_aligned.bam shortReads_chimp1_aligned.sorted
samtools index shortReads_chimp1_aligned.sorted.bam
```
