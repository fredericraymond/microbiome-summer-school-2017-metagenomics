# microbiome-summer-school-2017-metagenomics

This is the tutorial for the assembly-based metagenomic course of the 2017 Microbiome Summer School.

## Software requirements

For this tutorial, you will need the following software :

* Ray Meta : https://github.com/zorino/ray
* Ray Platform (required to install Ray) : https://github.com/sebhtml/RayPlatform
* Ray Meta Tools : https://github.com/fredericraymond/RayMetaTools
* MPI : https://www.open-mpi.org/software/ompi/v2.1/

The tutorial will be performed on the following dataset :

TBD

## Introduction

Ray Meta is one of the many softwares that allow the assembly of metagenomes. As explained in the plenary session, each assembly software has its advantages and inconveniences. The main advantage of Ray Meta is that it integrates assembly with taxonomical profiling. It also allow to identify the putative taxonomical origin of contigs. These are tools that are very useful when working with metagenomic data. The scalability of Ray makes it supercomputer friendly and thus useful to work on large metagenomic studies with samples that have hundreds of millions of reads. As a drawback, Ray Meta may output smaller assemblies that other softwares such as MegaHit or MetaSpades. During your metagenomic work, it is important to consider all the available tools in order to make the best of your data. In this tutorial, I will show how to use Ray Meta and interpret its output in order to better understand important concepts relating to assembly-based metagenomics.

## Step 1 - Assembling a microbiome using Ray Meta

In the first part of this tutorial, we will assemble a mixed bacterial culture and evaluate the quality of the assembly.

For the purpose of this tutorial, we will use a rarefied dataset that will not provide a complete assembly and then compare the results with a sample with sufficient sequencing depth.

### Assemble and profile sample

The command to run the assembly and profiling using Ray Meta is as follow :


```
mpiexec -n 2 Ray \
 -o \
 Sample_RVH-2106-Ray-2017-06-18 \
 -k \
 21 \
 -p \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R1_001.fastq.gz \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R2_001.fastq.gz \
 -search \
 genomes \
 -with-taxonomy \
 Genome-to-Taxon.tsv \
 TreeOfLife-Edges.tsv \
 Taxon-Names.tsv
```

Please go to the ~/SummerSchoolMicrobiome/UsingRayMeta/ folder and run the preceding command. 


```
cd ~/SummerSchoolMicrobiome/UsingRayMeta/

mpiexec -n 2 Ray \
 -o \
 Sample_RVH-2106-Ray-2017-06-18 \
 -k \
 21 \
 -p \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R1_001.fastq.gz \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R2_001.fastq.gz \
 -search \
 genomes \
 -with-taxonomy \
 Genome-to-Taxon.tsv \
 TreeOfLife-Edges.tsv \
 Taxon-Names.tsv
```

This will take around 45 minutes, so we will continue the tutorial below. Now that it is running, We will examine this command line-by-line.

```
mpiexec -n 2 Ray \
```

Mpiexec is the paralellization software used by Ray. The "-n" allows to indicate how many processors should be used for the analysis. For this tutorial, we will be using 2 cores. In a real project, the number of processors will depend on the amount of sequence to be assembled and on the size of the profiling datasets (which we will explain soon). To assemble a single genome, we often use a value of 32 while for complex microbiomes we can use up to 128 and 256 processors.

```
 -o \
 Sample_RVH-2106-Ray-2017-06-06 \
```

The "-o" is the output directory.


```
 -k \
 31 \
```
 
The "-k" is the k-mer length used for the analysis. In a majority of cases a k-mer length of 31 will provide optimal assembly. However, depending on the sequencing depth and on the complexity of the organism or community, a different value for k could be optimal. In doubt, it can be useful to compare assemblies using different values of k.
 
```
 -p \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R1_001.fastq.gz \
 Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R2_001.fastq.gz \
```

The "-p" tells Ray to use paired-end reads. Therefore, we give it to files, one for read 1 and one for read to. Alternatively, "-s" can be used for single reads.

```
 -search \
 genomes \
```

The "-search" function allows to quantify the k-mer content based on a reference database. In this example, the genome directory contains three genomes :

* Clostridium_difficile_630.fasta
* Escherichia_coli_K_12_substr__DH10B.fasta
* Leuconostoc_lactis_KCTC_3528_uid68683.fasta

In a real analysis, this directory can contain thousands of bacterial genomes. Several directories can be given to the software by repeating the "-search" argument. For example, if working on the human microbiome, the following sequence file could be included in addition to bacterial genomes : "-search humangenome". This is important as k-mers from the host can sometimes cause false positive signals for certain bacteria.

```
 -with-taxonomy \
 Genome-to-Taxon.tsv \
 TreeOfLife-Edges.tsv \
 Taxon-Names.tsv
```
This final argument tells Ray which taxonomy to use.

* *Genome-to-Taxon.tsv* : makes the link between genome sequences and their taxonomical classification. For this tutorial, we use a simplified version for the Genome-to-Taxon.tsv file since the complete version can be several Gb. 
* *TreeOfLife-Edges.tsv* : Allows to reconstitute the taxonomical tree.
* *Taxon-Names.tsv* : Gives the names of the nodes of the taxonomical tree (the taxa names).

When adding new genomes to the taxonomy, one must make sure to include the matching information in these three files.


## Step 2 - Understanding the Ray output directory

How is our assembly doing? To see at which step the assembly is currently, we can consult the file *ElapsedTime.txt*.

```
cd ~/SummerSchoolMicrobiome/UsingRayMeta/Sample_RVH-2106-Ray-2017-06-18
more ElapsedTime.txt
```

Precomputed results for the assembly of the sample is available in the following directory :  *~/SummerSchoolMicrobiome/UsingRayMeta/Sample_RVH-2106-Ray-2017-06-06*.

We will look at the final version of the precomputed assembly. Note that this assembly was done using 4 processors, so it took approximately half the time of our current assembly.

```
cd ~/SummerSchoolMicrobiome/UsingRayMeta/Sample_RVH-2106-Ray-2017-06-06
more ElapsedTime.txt

  Network testing	2017-06-06T15:17:42	0 seconds	0 seconds
  Counting sequences to assemble	2017-06-06T15:17:44	2 seconds	2 seconds
  Sequence loading	2017-06-06T15:17:52	8 seconds	10 seconds
  K-mer counting	2017-06-06T15:18:18	26 seconds	36 seconds
  Coverage distribution analysis	2017-06-06T15:18:21	3 seconds	39 seconds
  Graph construction	2017-06-06T15:18:53	32 seconds	1 minutes, 11 seconds
  Null edge purging	2017-06-06T15:19:35	42 seconds	1 minutes, 53 seconds
  Selection of optimal read markers	2017-06-06T15:20:20	45 seconds	2 minutes, 38 seconds
  Detection of assembly seeds	2017-06-06T15:21:57	1 minutes, 37 seconds	4 minutes, 15 seconds
  Estimation of outer distances for paired reads	2017-06-06T15:22:02	5 seconds	4 minutes, 20 seconds
  Bidirectional extension of seeds	2017-06-06T15:23:20	1 minutes, 18 seconds	5 minutes, 38 seconds
  Merging of redundant paths	2017-06-06T15:24:10	50 seconds	6 minutes, 28 seconds
  Generation of contigs	2017-06-06T15:24:11	1 seconds	6 minutes, 29 seconds
  Scaffolding of contigs	2017-06-06T15:25:19	1 minutes, 8 seconds	7 minutes, 37 seconds
  Counting sequences to search	2017-06-06T15:25:19	0 seconds	7 minutes, 37 seconds
  Graph coloring	2017-06-06T15:25:28	9 seconds	7 minutes, 46 seconds
  Counting contig biological abundances	2017-06-06T15:25:31	3 seconds	7 minutes, 49 seconds
  Counting sequence biological abundances	2017-06-06T15:25:39	8 seconds	7 minutes, 57 seconds
  Loading taxons	2017-06-06T15:25:42	3 seconds	8 minutes, 0 seconds
  Loading tree	2017-06-06T15:25:51	9 seconds	8 minutes, 9 seconds
  Processing gene ontologies	2017-06-06T15:25:57	6 seconds	8 minutes, 15 seconds
  Computing neighbourhoods	2017-06-06T15:25:57	0 seconds	8 minutes, 15 seconds

```

Now, is this assembly any good? To investigate, we will look at the *OutputNumbers.txt* file.

```
more OutputNumbers.txt

Contigs >= 100 nt
 Number: 9649
 Total length: 1862984
 Average: 193
 N50: 196
 Median: 148
 Largest: 6958
Contigs >= 500 nt
 Number: 305
 Total length: 253235
 Average: 830
 N50: 764
 Median: 644
 Largest: 6958
Scaffolds >= 100 nt
 Number: 9557
 Total length: 1880491
 Average: 196
 N50: 197
 Median: 147
 Largest: 6958
Scaffolds >= 500 nt
 Number: 251
 Total length: 288026
 Average: 1147
 N50: 1449
 Median: 707
 Largest: 6958
```

<details> 
  <summary>Q1: What can you say about his assembly? Is it good? Do we have good coverage of the expected genome? </summary>
   The total length of this assembly is low, less than 2 million including all contigs. When looking only at contigs longer than 500 nucleotide, the assembly is only 250 kb. Here, we were assembling a Clostridium difficile genome. Thus, we were expecting a genome length around 4 Mb.
</details>

To troubleshoot these results,








What's important and what's cool

## Step 3 - Profiling results

Yep, make a big table

## Step 4 - Contig identification

Bin it baby
