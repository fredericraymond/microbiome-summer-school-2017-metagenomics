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


How do you interpret these statistics?


<details> 
  <summary>Q1: What can you say about his assembly? Is it good? Do we have good coverage of the expected genome? </summary>
   The total length of this assembly is low, less than 2 million including all contigs. When looking only at contigs longer than 500 nucleotide, the assembly is only 250 kb.  The number of contigs is very high for the length of the assembly. Here, we were assembling a Clostridium difficile genome. Thus, we were expecting a genome length around 4 Mb. To troubleshoot these results, we will look at the sequencing statistics.
</details>
<br>
<br>
The number of reads used for assembly can be found in the file *NumberOfSequences.txt*


```
more NumberOfSequences.txt

  Files: 2

FileNumber: 0
	FilePath: Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R1_001.fastq.gz
 	NumberOfSequences: 1341231
	FirstSequence: 0
	LastSequence: 1341230

FileNumber: 1
	FilePath: Sample_RVH-2106/RVH-2106_GTAGAGGA-TAGATCGC_L003_R2_001.fastq.gz
 	NumberOfSequences: 1341231
	FirstSequence: 1341231
	LastSequence: 2682461


Summary
	NumberOfSequences: 2682462
	FirstSequence: 0
	LastSequence: 2682461

```

The number of reads in the analysis is 2,682,462. The number of nucleotides sequenced is 37 * 2,682,462 = 270,928,662. If we divide it by the expected genome length, we get 270,928,662 / 4e6 = 68X coverage, which should be sufficient for analysis. However, it is obviously not the case.

Further troubleshooting can be done by looking at the *CoverageDistributionAnalysis.txt* file 

```
more CoverageDistributionAnalysis.txt

k-mer length:   21
Number of k-mers in the distributed de Bruijn graph: 10526772
Lowest coverage observed:       2
MinimumCoverage:        2
PeakCoverage:   2
RepeatCoverage: 2
Number of k-mers with at least MinimumCoverage: 10526772 k-mers
Percentage of vertices with coverage 2: 12.1506 %
DistributionFile: Sample_RVH-2106-Ray-2017-06-06/CoverageDistribution.txt
```

The information of interest is the number of k-mers with at least MinimumCoverage: 10526772 k-mers. This estimates the size of the de bruijn graph that attains the minimum coverage, in this case a depth of 2. If we divide this number by 2 (to account for the two strands of DNA), we can estimate the actual size of the genome or metagenome. In this case, the genome size would be approximately 5 Mb, which is slightly bigger than a C. difficile genome.

We will look into to profiling results of Ray to investigate the taxonomical content of the sample. This information is found in the BiologicalAbundance folder.

```
cd BiologicalAbundances
ls

0.Profile.genomes.tsv              0.Profile.TaxonomyRank=no rank.tsv  0.Profile.TaxonomyRank=superkingdom.tsv  _Frequencies
0.Profile.TaxonomyRank=class.tsv   0.Profile.TaxonomyRank=order.tsv    _Coloring                                _GeneOntology
0.Profile.TaxonomyRank=family.tsv  0.Profile.TaxonomyRank=phylum.tsv   _DeNovoAssembly                          genomes
0.Profile.TaxonomyRank=genus.tsv   0.Profile.TaxonomyRank=species.tsv  _Directories.tsv                         _Taxonomy

```

The .tsv files provide the proportion of samples that could be identified as taxonomic units (in k-mers).

<details> 
  <summary>Q2: Take a few minutes to look at the .tsv files and try to understand what is happening with this assembly.</summary>
   In *0.Profile.TaxonomyRank=species.tsv*, we can see that 35% of the sample is *Leuconostoc lactis* while the expected *C. difficile* is only 65% of the sample. The *E. coli* signal is background noise. Yep, this sample is not pure. It is a mixed culture.
</details>
<br>
<br>
To simplify the interpretation when working with mixed population, we can look at the taxa that are at more than 0.01% of a sample.In this example, it is not that useful, but as you can see it removes the background signal from *E. coli*.

```
 awk -F "\t" '$4>0.0001' 0.Profile.TaxonomyRank=species.tsv

  #TaxonIdentifier        TaxonName       TaxonRank       TaxonProportion
  1246    Leuconostoc lactis      species 0.352636
  1496    Peptoclostridium difficile      species 0.647357
```

Now that we have troubleshooted our sample, we will use an assembly of the same sample with more reads to look at a proper assembly of these genomes.

Now, we will move to another directory. In this second assembly, take a look at the same files we have consulted in this section to answer the following questions :

```
cd ~/SummerSchoolMicrobiome/UsingRayMeta/Sample_RVH-2106-Ray-2017-06-18
```

How many reads did we use?
What is the assembly size? Broken in how many contigs?
What is is the size of the graph? How does it relate with the size of the assembly and the actual size of a *C. difficile* genome?
How is the species proportion compared to the lower sequencing depth analysis?

## Step 3 - Contig binning

From this complete assembly, we will bin the contigs to separate the *Leuconostoc* from the *Clostridium*.

Some data on contig identification is already available in the files, although we will need to do some processing to make this information more valuable.

The first file to look into is in the *BiologicalAbundances/genomes/ContigIdentification.tsv*.

```
more BiologicalAbundances/genomes/ContigIdentification.tsv

#Contig name    K-mer length    Contig length in k-mers Contig strand   Category        Sequence number Sequence name   Sequence length in k-mers       Matches in contig       Contig length ratio     Sequence length ratio
contig-573000003        21      235     F       Leuconostoc_lactis_KCTC_3528_uid68683   0       gi|339303864|ref|NZ_AEOR01000001.1| Leuconostoc lactis KCTC 3528        882     142     0.604255        0.160998
contig-743000002        21      207     F       Leuconostoc_lactis_KCTC_3528_uid68683   0       gi|339303864|ref|NZ_AEOR01000001.1| Leuconostoc lactis KCTC 3528        882     157     0.758454        0.178005
contig-494000000        21      183     R       Leuconostoc_lactis_KCTC_3528_uid68683   0       gi|339303864|ref|NZ_AEOR01000001.1| Leuconostoc lactis KCTC 3528        882     82      0.448087        0.0929705
contig-1276000002       21      346     R       Leuconostoc_lactis_KCTC_3528_uid68683   0       gi|339303864|ref|NZ_AEOR01000001.1| Leuconostoc lactis KCTC 3528        882     224     0.647399        0.253968
contig-778000001        21      124     F       Leuconostoc_lactis_KCTC_3528_uid68683   2       gi|339303866|ref|NZ_AEOR01000003.1| Leuconostoc lactis KCTC 3528        1060    119     0.959677        0.112264
contig-1972000000       21      146     F       Leuconostoc_lactis_KCTC_3528_uid68683   2       gi|339303866|ref|NZ_AEOR01000003.1| Leuconostoc lactis KCTC 3528        1060    106     0.726027        0.1
contig-1065000002       21      123     R       Leuconostoc_lactis_KCTC_3528_uid68683   2       gi|339303866|ref|NZ_AEOR01000003.1| Leuconostoc lactis KCTC 3528        1060    123     1       0.116038
contig-327000001        21      319     R       Leuconostoc_lactis_KCTC_3528_uid68683   4       gi|339303868|ref|NZ_AEOR01000005.1| Leuconostoc lactis KCTC 3528        654     319     1       0.487768
contig-553000000        21      732     R       Leuconostoc_lactis_KCTC_3528_uid68683   4       gi|339303868|ref|NZ_AEOR01000005.1| Leuconostoc lactis KCTC 3528        654     169     0.230874        0.25841
```

For each contig, we get several information, notably the lenght of the contig, the reference sequence matching the contig and the proportion of the contig matched by the reference sequence. The problem (or the advantage?) with this file is that it includes several entry for each contig, providing its similarity with several sequences from the reference database.  For a small community like this one, we could analyse this file. However, for more complex metagenome, we need help!


We will not be running the following steps since it take a while. It has already been computed and available in the directory. 

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
Contig-name     Contig_length_in_kmers  Contig_mode_kmer_depth  Total_colored_kmer      LCA_taxon_id    LCA_name        LCA_rank        LCA_score       phylum  phylum_score    class   class_score     order   order_score     family  family_score       genus   genus_score     species species_score
contig-31000024 322     14      287     1358    Lactococcus lactis      species 0.483870967742  Firmicutes      0.483870967742  Bacilli 0.483870967742  Lactobacillales 0.483870967742  Streptococcaceae        0.483870967742  Lactococcus        0.483870967742  Lactococcus lactis      0.483870967742
contig-31000025 242     6       0       0       No name Unknown rank    0               0               0               0               0               0               0
contig-31000026 148     2       4       1578    Lactobacillus   genus   1.0     Firmicutes      1.0     Bacilli 1.0     Lactobacillales 1.0     Lactobacillaceae        1.0     Lactobacillus 1.0     Lactobacillus zeae      0.5
contig-31000027 146     6       127     1358    Lactococcus lactis      species 1.0     Firmicutes      1.0     Bacilli 1.0     Lactobacillales 1.0     Streptococcaceae        1.0     Lactococcus     1.0     Lactococcus lactis      1.0
contig-31000020 353     11      16      unknown No name Unknown rank    0               -1              -1              -1              -1              -1              -1
contig-31000021 225     7       0       0       No name Unknown rank    0               0               0               0               0               0               0
contig-31000022 162     2       80      151534  Lactococcus phage bIL311        species 1.0             -1              -1      Caudovirales    1.0     Siphoviridae    1.0             -1      Lactococcus phage bIL311        1.0
contig-31000023 165     5       0       0       No name Unknown rank    0               0               0               0               0               0               0
contig-31000028 159     3       159     186826  Lactobacillales order   0.521869158879  Firmicutes      0.780934579439  Bacilli 0.521869158879  Lactobacillales 0.521869158879  Lactobacillaceae        0.36953271028   Lactobacillus   0.338878504673     Clostridium butyricum   0.135514018692
```

In this file, we put all the beauty of Ray Meta. For each contig, we get the following information :

* Contig-name
* Contig_length_in_kmers  
* Contig_mode_kmer_depth  - Mode coverage depth of the k-mers in the contig.
* Total_colored_kmer      - Number of k-mers in the contig that have been colored by the reference database.
* LCA_taxon_id            - Best hit taxonomical association NCBI taxID.
* LCA_name                - Best hit taxonomical association taxa name.
* LCA_rank                - Best hit taxonomical association; rank of the best hit.
* LCA_score               - Best hit taxonomical association score

The following colums indicate the best classification at each taxonomical rank and its score. A score > 0.20 indicate a high probability that the taxonomical association is good.

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
* species species_score

The score, also known as the pl-value, is computed using the following formula :

###

We can now bin contigs and compute the size of each bin that is associated to a species. We will work with R.

Let's open R.

```
R
```

Load the contig identification file.

```
data <- read.table("Sample_503619-Ray-2016-12-10_original-nemesis.tsv", header=1, row.names=1, sep="\t")
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

How many families have been identified in our file ?

```
summary(data$family)
```

We select the lines for which species are identified with a score greater than 0.2.

####



Analysis of laughing-nemesis results.

Explanation of LCA.

Working with the file to sort it.




