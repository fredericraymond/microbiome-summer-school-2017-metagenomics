# Microbiome Summer School 2017 - Metagenomics using Ray Meta

This is the tutorial for the assembly-based metagenomic course of the 2017 Microbiome Summer School.

## Software requirements

For this tutorial, you will need the following software :

* Ray Meta : https://github.com/zorino/ray
* Ray Platform (required to install Ray) : https://github.com/sebhtml/RayPlatform
* MPI : https://www.open-mpi.org/software/ompi/v2.1/
* R : https://www.r-project.org/
* Laughing-nemesis : https://github.com/plpla/laughing-nemesis.git
* Python

The tutorial will be performed on the following dataset :

TBD

## Introduction

Ray Meta is one of the many softwares that allow the assembly of metagenomes. As explained in the plenary session, each assembly software has its advantages and inconveniences. The main advantage of Ray Meta is that it integrates assembly with taxonomical profiling. It also allows to identify the putative taxonomical origin of contigs. These are tools that are very useful when working with metagenomic data. The scalability of Ray makes it supercomputer friendly and thus useful to work on large metagenomic studies with samples that have hundreds of millions of reads. As a drawback, current versions Ray Meta may output smaller assemblies that other softwares such as MegaHit or MetaSpades. During your metagenomic work, it is important to consider all the available tools in order to make the best of your data. In this tutorial, I will show how to use Ray Meta and interpret its output in order to better understand important concepts relating to assembly-based metagenomics.

## Step 0 - Download data

For this tutorial, we will use a syntetic dataset inspired by a genome we once sequenced for a collaborator. So, in theory, we are sequencing a genome from *Clostridium difficile*...

```
wget xxx
```

## Step 1 - Metagenome assembly using Ray Meta

In the first part of this tutorial, we will assemble a mixed bacterial culture and evaluate the quality of the assembly.

For the purpose of this tutorial, we will use a rarefied dataset that will not provide a complete assembly and then compare the results with a sample with sufficient sequencing depth.

### Assemble and profile sample

The command to run the assembly and profiling using Ray Meta is as follow. Please go to the ~/UsingRayMeta/ folder and run the preceding command. 


```
cd ~/UsingRayMeta/

mpiexec -n 2 Ray \
 -o \
 Sample-Small \
 -k \
 31 \
 -p \
 reads/MetaSim_Cdiff8_Ecoli2_50000_fir.fastq \
 reads/MetaSim_Cdiff8_Ecoli2_50000_sec.fastq \
 -search \
 genomes \
 -with-taxonomy \
 Genome-to-Taxon.tsv \
 TreeOfLife-Edges.tsv \
 Taxon-Names.tsv
```

This will take around 45 minutes, so we will continue the tutorial below. Now that it is running, we will examine this command line-by-line.

```
mpiexec -n 2 Ray \
```

Mpiexec is the paralellization software used by Ray. The "-n" allows to indicate how many cores should be used for the analysis. For this tutorial, we will be using 2 cores. In a real project, the number of processors will depend on the amount of sequence to be assembled and on the size of the profiling datasets (which we will explain soon). To assemble a single genome, we often use a value of 32, while for complex microbiomes we can use up to 128 or 256 processors.

```
 -o \
 Sample-Small \
```

The "-o" is the output directory.


```
 -k \
 31 \
```
 
The "-k" is the k-mer length used for the analysis. In a majority of cases a k-mer length of 31 will provide optimal assembly. However, depending on the sequencing depth and on the complexity of the organism or community, a different value for k could be optimal. In doubt, it can be useful to compare assemblies using different values of k. In our lab, it is standard procedure to test 21, 31, 51 and 71 as values of k.
 
```
 -p \
 reads/MetaSim_Cdiff8_Ecoli2_50000_fir.fastq \
 reads/MetaSim_Cdiff8_Ecoli2_50000_sec.fastq \
```

The "-p" tells Ray to use paired-end reads. Therefore, we give it two files, one for read 1 and one for read 2. Alternatively, "-s" can be used for single reads. Several entries of -p or -s can be added to the arguments

```
 -search \
 genomes \
```

The "-search" function allows to quantify the k-mer content based on a reference database. In this example, the genome directory contains three genomes :

* Clostridium_difficile_630.fasta
* Escherichia_coli_K_12_substr__DH10B.fasta
* Leuconostoc_lactis_KCTC_3528_uid68683.fasta

In a real analysis, this directory can contain thousands of bacterial genomes. Several directories can be given to the software by repeating the "-search" argument. For example, if working on the human microbiome, the following sequence file could be included in addition to bacterial genomes (could be : "-search GRCh38"). This is important as k-mers from the host can sometimes cause false positive signals for certain bacteria.

```
 -with-taxonomy \
 Genome-to-Taxon.tsv \
 TreeOfLife-Edges.tsv \
 Taxon-Names.tsv
```

This final argument tells Ray which taxonomy to use.

* *Genome-to-Taxon.tsv* : Makes the link between genome sequences and their taxonomical classification. For this tutorial, we use a simplified version of the Genome-to-Taxon.tsv file, since the complete version can be several Gb. 
* *TreeOfLife-Edges.tsv* : Allows to reconstitute the taxonomical tree.
* *Taxon-Names.tsv* : Gives the names of the nodes of the taxonomical tree (the taxa names).

When adding new genomes to the taxonomy, one must make sure to include the matching information in these three files.


## Step 2 - Understanding the Ray output directory

How is our assembly doing? To see at which step the assembly is currently, we can consult the file *ElapsedTime.txt*.

```
cd ~/UsingRayMeta/Sample-Small
more ElapsedTime.txt
```

Precomputed results for the assembly of this sample (but with a lot more reads!) is available in the following directory : *~/UsingRayMeta/Sample-Big*.

We will look at the final version of the precomputed assembly. Note that this assembly was done using 48 processors, so it took only 20 minutes to run.

```
cd ~/UsingRayMeta/Sample-Big
more ElapsedTime.txt

  #Step	Date	Elapsed time	Since Beginning
  Network testing	2017-06-15T14:47:08	1 seconds	1 seconds
  Counting sequences to assemble	2017-06-15T14:47:09	1 seconds	2 seconds
  Sequence loading	2017-06-15T14:47:16	7 seconds	9 seconds
  K-mer counting	2017-06-15T14:47:36	20 seconds	29 seconds
  Coverage distribution analysis	2017-06-15T14:47:43	7 seconds	36 seconds
  Graph construction	2017-06-15T14:48:09	26 seconds	1 minutes, 2 seconds
  Null edge purging	2017-06-15T14:48:49	40 seconds	1 minutes, 42 seconds
  Selection of optimal read markers	2017-06-15T14:49:24	35 seconds	2 minutes, 17 seconds
  Detection of assembly seeds	2017-06-15T14:51:39	2 minutes, 15 seconds	4 minutes, 32 seconds
  Estimation of outer distances for paired reads	2017-06-15T14:51:52	13 seconds	4 minutes, 45 seconds
  Bidirectional extension of seeds	2017-06-15T14:55:36	3 minutes, 44 seconds	8 minutes, 29 seconds
  Merging of redundant paths	2017-06-15T15:03:44	8 minutes, 8 seconds	16 minutes, 37 seconds
  Generation of contigs	2017-06-15T15:03:53	9 seconds	16 minutes, 46 seconds
  Scaffolding of contigs	2017-06-15T15:05:33	1 minutes, 40 seconds	18 minutes, 26 seconds
  Counting sequences to search	2017-06-15T15:05:34	1 seconds	18 minutes, 27 seconds
  Graph coloring	2017-06-15T15:06:01	27 seconds	18 minutes, 54 seconds
  Counting contig biological abundances	2017-06-15T15:06:13	12 seconds	19 minutes, 6 seconds
  Counting sequence biological abundances	2017-06-15T15:06:55	42 seconds	19 minutes, 48 seconds
  Loading taxons	2017-06-15T15:07:02	7 seconds	19 minutes, 55 seconds
  Loading tree	2017-06-15T15:07:22	20 seconds	20 minutes, 15 seconds
  Processing gene ontologies	2017-06-15T15:07:34	12 seconds	20 minutes, 27 seconds
  Computing neighbourhoods	2017-06-15T15:07:34	0 seconds	20 minutes, 27 seconds

```

Now, is this assembly any good? To investigate, we will look at the *OutputNumbers.txt* file.

```
more OutputNumbers.txt

Contigs >= 100 nt
 Number: 13633
 Total length: 8137962
 Average: 596
 N50: 10391
 Median: 234
 Largest: 265237
Contigs >= 500 nt
 Number: 1531
 Total length: 5224457
 Average: 3412
 N50: 85464
 Median: 632
 Largest: 265237
Scaffolds >= 100 nt
 Number: 13619
 Total length: 8139296
 Average: 597
 N50: 15499
 Median: 234
 Largest: 265237
Scaffolds >= 500 nt
 Number: 1517
 Total length: 5225791
 Average: 3444
 N50: 109710
 Median: 631
 Largest: 265237

```


How do you interpret these statistics?


<details> 
  <summary>Q1: What can you say about his assembly? Is it good? Do we have good coverage of the expected C. difficile genome? </summary>
   The total length of this assembly is higher than what is expecter for a C. difficile genome. When looking only at contigs longer than 500 nucleotide, the assembly is 5,224,457 nt.  The number of contigs is very high for the length of the assembly. Here, we were assembling a Clostridium difficile genome. Thus, we were expecting a genome length around 4 Mb. To troubleshoot these results, we will look at the sequencing statistics.
</details>
<br>
<br>
The number of reads used for assembly can be found in the file *NumberOfSequences.txt*


```
more NumberOfSequences.txt

  Files: 2

FileNumber: 0
	FilePath: MetaSim_Cdiff8_Ecoli2_1000000_fir.fastq
 	NumberOfSequences: 1000000
	FirstSequence: 0
	LastSequence: 999999

FileNumber: 1
	FilePath: MetaSim_Cdiff8_Ecoli2_1000000_sec.fastq
 	NumberOfSequences: 1000000
	FirstSequence: 1000000
	LastSequence: 1999999


Summary
	NumberOfSequences: 2000000
	FirstSequence: 0
	LastSequence: 1999999

```

The number of reads in the analysis is 2,000,000. The number of nucleotides sequenced is 101 * 2,000,000 = 202,000,000. If we divide it by the expected genome length, we get 202,000,000 / 4e6 = 50X coverage, which is indeed sufficient for assembly.

Further troubleshooting can be done by looking at the *CoverageDistributionAnalysis.txt* file 

```
more CoverageDistributionAnalysis.txt

k-mer length:	31
Number of k-mers in the distributed de Bruijn graph: 18616594
Lowest coverage observed:	2
MinimumCoverage:	10
PeakCoverage:	17
RepeatCoverage:	24
Number of k-mers with at least MinimumCoverage:	8459298 k-mers
Percentage of vertices with coverage 2:	15.5668 %
DistributionFile: Sample_SS-Big/CoverageDistribution.txt

```

The information of interest is the number of k-mers with at least MinimumCoverage: 8,459,298 k-mers. This estimates the size of the de bruijn graph that attains the minimum coverage, in this case a depth of 10. If we divide this number by 2 (to account for the two strands of DNA), we can estimate the actual size of the genome or metagenome. In this case, the genome size would be approximately 4,2 Mb, which is what we would expect as a genome size. However, the number of k-mers in the graph is almost twice that size : 18,616,594. This may indicate something is not behaving as expected...

We will look into to profiling results of Ray to investigate the taxonomical content of the sample. This information is found in the BiologicalAbundance folder.

```
cd BiologicalAbundances
ls

0.Profile.genomes.tsv              0.Profile.TaxonomyRank=no rank.tsv  0.Profile.TaxonomyRank=superkingdom.tsv  _Directories.tsv  _Taxonomy
0.Profile.TaxonomyRank=class.tsv   0.Profile.TaxonomyRank=order.tsv    20170616_Nemesis.tsv                     _Frequencies
0.Profile.TaxonomyRank=family.tsv  0.Profile.TaxonomyRank=phylum.tsv   _Coloring                                _GeneOntology
0.Profile.TaxonomyRank=genus.tsv   0.Profile.TaxonomyRank=species.tsv  _DeNovoAssembly                          genomes

```

The .tsv files provide the proportion of samples that could be identified as taxonomic units (in k-mers).

<details> 
  <summary>Q2: Take a few minutes to look at the .tsv files and try to understand what is happening with this assembly.</summary>
   In *0.Profile.TaxonomyRank=species.tsv*, we can see that 20% of the sample is *Leuconostoc lactis* while the expected *C. difficile* is only 20% of the sample. The *E. coli* signal is background noise (< 10e-6). Yep, this sample is not pure. It is a mixed culture.
</details>
<br>
<br>
To simplify the interpretation when working with mixed population, we can look at the taxa that are at more than 0.01% of a sample. In this example, it is not that useful, but as you can see it removes the background signal from *E. coli*.

```
 awk -F "\t" '$4>0.0001' 0.Profile.TaxonomyRank=species.tsv

  #TaxonIdentifier        TaxonName       TaxonRank       TaxonProportion
  562     Escherichia coli        species 0.209279
  1496    Peptoclostridium difficile      species 0.79072

```

Now, we will move to another directory, with more complex metagenomic sample. In this second assembly, take a look at the same files we have consulted in this section to answer the following questions :

```
cd ~/UsingRayMeta/Sample_P4J7-Assembly
```

How many reads did we use?
What is the assembly size? Broken in how many contigs?
What is is the size of the graph? How does it relate with the size of the assembly?
How is the species proportion compared to the lower sequencing depth analysis?

## Step 3 - Contig binning

From this complete assembly, we will bin the contigs to separate the *Escherichia* from the *Clostridium*.

Some data on contig identification is already available in the files, although we will need to do some processing to make this information more valuable.

The first file to look into is in the *BiologicalAbundances/genomes/ContigIdentification.tsv*.

```
cd ~/UsingRayMeta/Sample-Big

more BiologicalAbundances/genomes/ContigIdentification.tsv

#Contig name	K-mer length	Contig length in k-mers	Contig strand	Category	Sequence number	Sequence name	Sequence length in k-mers	Matches in contig	Contig length ratio	Sequence length ratio
contig-0	31	37550	F	Clostridium_difficile_630	1	gi|126697566|ref|NC_009089.1| Clostridium difficile 630, complet	4290222	103	0.00274301	2.40E-05
contig-0	31	37550	R	Clostridium_difficile_630	1	gi|126697566|ref|NC_009089.1| Clostridium difficile 630, complet	4290222	37401	0.996032	0.00871773
contig-1	31	1601	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	1574	0.983136	0.000335886
contig-1	31	1601	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	27	0.0168645	5.76E-06
contig-10	31	1401	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	1401	1	0.000298969
contig-1000000	31	82506	F	Clostridium_difficile_630	1	gi|126697566|ref|NC_009089.1| Clostridium difficile 630, complet	4290222	82447	0.999285	0.0192174
contig-10000000	31	717	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	717	1	0.000153005
contig-100000000	31	276	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	276	1	5.89E-05
contig-100000001	31	234	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	234	1	4.99E-05
contig-100000002	31	254	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	254	1	5.42E-05
contig-100000003	31	256	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	256	1	5.46E-05
contig-100000004	31	265	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	265	1	5.66E-05
contig-100000005	31	248	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	248	1	5.29E-05
contig-100000006	31	273	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	273	1	5.83E-05
contig-100000007	31	267	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	267	1	5.70E-05
contig-100000008	31	261	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	261	1	5.57E-05
contig-100000009	31	269	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	269	1	5.74E-05
contig-10000001	31	687	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	687	1	0.000146604
contig-100000010	31	239	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	239	1	5.10E-05
contig-100000011	31	259	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	259	1	5.53E-05
contig-100000012	31	1737	R	Leuconostoc_lactis_KCTC_3528_uid68683	1149	gi|339305013|ref|NZ_AEOR01001150.1| Leuconostoc lactis KCTC 3528	2075	16	0.00921128	0.00771084
contig-100000012	31	1737	F	Clostridium_difficile_630	1	gi|126697566|ref|NC_009089.1| Clostridium difficile 630, complet	4290222	35	0.0201497	8.16E-06
contig-100000012	31	1737	F	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	1620	0.932642	0.000345703
contig-100000012	31	1737	R	Escherichia_coli_K_12_substr__DH10B	0	gi|170079663|ref|NC_010473.1| Escherichia coli str. K-12 substr.	4686107	117	0.0673575	2.50E-05

```

For each contig, we get several information, notably the lenght of the contig, the reference sequence matching the contig and the proportion of the contig matched by the reference sequence. The problem (or the advantage?) with this file is that it includes several entry for each contig (see contig-100000012), providing its similarity with several sequences from the reference database.  For a small community like this one, we could analyse this file. However, for more complex metagenome, we need help!

We will not be running the following steps since it take some time. It has already been computed and available in the directory. 

Enter the following commands to download the laughing-nemesis last common ancestor tool for Ray Meta.

```
cd ~
mkdir software
cd software
git clone https://github.com/plpla/laughing-nemesis.git
```

We would then use tho following command line to run it on our assembly :

```
python ~/software/laughing-nemesis/FindLastCommonAncester.py lca -d ~/SummerSchoolMicrobiome/Sample_RVH-2106-Ray-2017-06-06/BiologicalAbundances -t ~/SummerSchoolMicrobiome/UsingRayMeta/TreeOfLife-Edges.tsv
```

Now, we look at the precomputed results.

```
more BiologicalAbundances/20170616_Nemesis.tsv

Contig-name	Contig_length_in_kmers	Contig_mode_kmer_depth	Total_colored_kmer	LCA_taxon_id	LCA_name	LCA_rank	LCA_score	phylum	phylum_score	class	class_score	order	order_score	family	family_score	genus	genus_score	species	species_score
contig-187000039	149	3	149	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-293000047	95	2	95	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-293000046	106	5	106	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-203000019	138	2	138	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-203000018	129	5	129	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-280000038	101	5	101	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-228000026	122	2	122	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-276000038	103	4	103	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-199000039	133	2	133	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-280000032	99	3	99	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-280000035	97	2	97	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-199000038	178	8	178	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-247000024	105	2	105	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-247000025	98	4	98	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-247000026	113	6	113	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1
contig-247000027	118	4	118	562	Escherichia coli	species	1	Proteobacteria	1	Gammaproteobacteria	1	Enterobacteriales	1	Enterobacteriaceae	1	Escherichia	1	Escherichia coli	1

```

In this file, we gather all the beauty of Ray Meta. For each contig, we get the following information :

* Contig-name
* Contig_length_in_kmers  
* Contig_mode_kmer_depth  - Mode coverage depth of the k-mers in the contig.
* Total_colored_kmer      - Number of k-mers in the contig that have been colored by the reference database.
* LCA_taxon_id            - Best hit taxonomical association NCBI taxID.
* LCA_name                - Best hit taxonomical association taxa name.
* LCA_rank                - Best hit taxonomical association rank of the best hit.
* LCA_score               - Best hit taxonomical association score

The following colums indicate the best classification at each taxonomical rank and its score. A score > 0.01 indicate a high probability that the taxonomical association is good.

* phylum  
* phylum_score    
* class   
* class_score     
* order   
* order_score     
* family  
* family_score       
* genus   
* genus_score     
* species
* species_score

The score, also known as the pl-value, is computed using the following formula :
Matches_in_contig^2 / (Contig_length_in_k-mers * Colored_k-mers)

We can now bin contigs and compute the size of each bin that is associated to a species. We will work with R.

Let's open R.

```
R
```

Load the contig identification file.

```
data <- read.table("~/UsingRayMeta/Sample-Big/BiologicalAbundances/20170616_Nemesis.tsv", header=1, row.names=1, sep="\t")
data[1:5,]
```

At which rank is the LCA the best ?

```
summary(data$LCA_rank)
```

How many families have been identified in our file ?

```
summary(data$family)
```

We select the lines for which species are identified with a score greater than 0.01.

```
data.select <- data[which(data$species_score>=0.01),]

### List species that were detected and sort them based on the number of contigs

sort(summary(droplevels(data.select$species)))
```

We calculate the sum of the length of contigs identified as Clostridium difficile (it is name *Peptoclostridium difficile* in this taxonomy).

```
sum(data.select[which(data.select$species=="Peptoclostridium difficile"),"Contig_length_in_kmers"])
```

We can generate a table to do that with all species that had at least one contig.

```
genome.size <- c()
for(species in levels(droplevels(data.select$species))){
  genome.size <- c(genome.size, sum(data.select[which(data.select$species==species),"Contig_length_in_kmers"]))
}
names(genome.size) <- levels(droplevels(data.select$species))
sort(genome.size)
```

We will look at the coverage distribution of contigs overall. This will allow to quality control contigs binning and better interpret results.

```
summary(data[,2])
plot(density(data[,2]))
```

Now, we look at a specific species.

```
summary(data.select[which(data.select$species=="Peptoclostridium difficile"),"Contig_mode_kmer_depth"])
plot(density(data.select[which(data.select$species=="Peptoclostridium difficile"),"Contig_mode_kmer_depth"]))

summary(data.select[which(data.select$species=="Escherichia coli"),"Contig_mode_kmer_depth"])
plot(density(data.select[which(data.select$species=="Escherichia coli"),"Contig_mode_kmer_depth"]))

```

We can also overlap the curve for the two species that we know to be present in this sample.

```
plot(density(data.select[which(data.select$species=="Peptoclostridium difficile"),"Contig_mode_kmer_depth"]))
lines(density(data.select[which(data.select$species=="Escherichia coli"),"Contig_mode_kmer_depth"]), col="red")
```

Finally, we may want generate a list of contigs that probably originate from SPECIES. We could with this list extract the contigs from the *Contigs.fasta* file in the assembly directory.

```
rownames(data.select[which(data.select$species=="Peptoclostridium difficile"),])
```

Want to try it with more complex data? Repeat the same procedure with THIS file.

```
data <- read.table("~/UsingRayMeta/Sample_P4J7-Assembly_nemesis.tsv", header=1, row.names=1, sep="\t")
data[1:5,]
```

Congratulation! You have completed the tutorial. Now, make good art!