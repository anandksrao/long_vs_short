# Comparing genome assemblies of a repetitive sequence using long and short reads
Genomes of many species including bacteria and mammals are rich in DNA sequences that are highly similar to other sequences in the genome. While many of these repetitive DNA segments seem to be non-functional others are thought to have unique functions and act as ‘selfish’ sequence modules. It is, therefore, of high importance to know the exact composition of these repeats.

Traditionally, these regions have been difficult to resolve due to the short read length of next generation sequencers. In de-novo assembly genomes may not be fully reconstructed if the reads are shorter than the repeats. In the alignment of reads to reference sequences repetitive segments of genomes are sometimes ignored to avoid false interpretation of the data. Thus, many recent genome assemblies have large structural gaps or misassembled regions in their sequences.

```
Several mechanisms have been shown to replicate and insert sequences in genomes. Repeats may be grouped into interspersed, tandem and nested repeats stretching from two to millions of nucleotides in length and copy number. In humans the two main classes of repeats are short tandem repeats (microsatellites) and long interspersed repeats (short/long interspersed nuclea elements; SINEs and LINEs).
```

It comes at no surprise then, that the negative impact of repeats on the accuracy of the assembly is largely avoided when read lengths are longer than the typical repeats.


### Aim
The aim of this tutorial is to demonstrate the effect of read length on assemblies of repetitive genome sequences. We will be comparing de-novo assemblies of two different sets of simulated reads: a set of short reads (~30 bases) and a set of long reads (~15 kb) and will calculate the proportion of the genome that was recovered by each.


### Further reading
- [Repetitive DNA and next-generation sequencing: computational challenges and solutions](http://www.nature.com/nrg/journal/v13/n1/full/nrg3117.html)
- [Reconstructing complex regions of genomes using long-read sequencing technology](http://genome.cshlp.org/content/24/4/688.full)

