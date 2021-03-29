---
layout: post
title: Siderastrea transcriptome mapping
date: '2021-01-25'
categories: Protocols
tags: [RNASeq, Bioinformatics]
---

# Coral transcriptome data processing
## ***Siderastrea***

### Analyzed by Veronica Radice, Barshis Lab, Old Dominion University

### 2-year transplant of ***Siderastrea siderea***
- Puerto Morelos, Mexico
- from low pH ojo (submarine discharge spring) and high pH control site to low pH ojo and control site
- Ana Martinez (PhD candidate) & Adina Paytan (PI) & Dan Barshis (PI)

### Transplant design (information from Ana)
- Corals were taken from 3 different origins: within the ojo (center), outside the ojo (control) and in the reef.
- Corals were then fixed in 3 transplant sites: within the ojo (Ojo LAJA Center and Ojo NORTE Center) and in a control site (Ojo Laja CONTROL).
- From each colony, 3 replicates (or individual cores) were taken

----------------------------------------------------------------------------------------------
## Data
- raw data (year 2 of transplant, fastq g-zipped) backed up on RC drive
- contains data for multiple species:  *Porites astreoides*, *Porites porites*, *Siderastrea siderea*
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2017_April/gslserver.qb3.berkeley.edu/170419_50SR_HS4K2A/Paytan

- other raw data from Paytan backed up on RC drive:
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/AdinaOA_2014_October_rawdata.tar.gz

### Files:
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/Siderastrea/


----------------------------------------------------------------------------------------------
## Notes from the literature - Symbiont genetics

*Siderastrea radians* symbiont associations
- Breviolum (ITS2 clade B)
  - Breviolum B5
    - B5 [https://link.springer.com/article/10.1007/s00227-020-03737-3](https://link.springer.com/article/10.1007/s00227-020-03737-3)
    - B5a (LaJeunesse 2002)
    - B5 [https://link.springer.com/article/10.1007%2Fs00338-006-0157-y](https://link.springer.com/article/10.1007%2Fs00338-006-0157-y)
    - B5 (Dr. Viridiana Avila-Magaña & Professor Mónica Medina, Pennsylvania State University) [http://medinalab.org/new/](http://medinalab.org/new/)
  - B1.I (70.31%) [https://link.springer.com/article/10.1007%2Fs00248-017-1096-6](https://link.springer.com/article/10.1007%2Fs00248-017-1096-6)
    - The four most abundant lineages in S. radians were B1.I, C1.I, B1.II, and C1.II (70.31, 13.41, 6.54, and 2.19%, respectively; Table S5; Fig. 3b)

*Siderastrea siderea* symbiont associations
- Breviolum
  - Breviolum psygmophilum (ITS2 type B2) 
    - [https://www.frontiersin.org/articles/10.3389/fmars.2018.00150/full](https://www.frontiersin.org/articles/10.3389/fmars.2018.00150/full)
- Cladocopium 
  - C1.I (74.39%) [https://link.springer.com/article/10.1007%2Fs00248-017-1096-6](https://link.springer.com/article/10.1007%2Fs00248-017-1096-6)
    - The four most abundant lineages in S. siderea were C1.I, C1.III, D1a, and B1.I (74.39, 12.94, 9.29, and 2.94%, respectively; Table S5; Fig. 3a)
  - C1 (LaJeunesse 2002)
  - C1 and C90 [https://peerj.com/articles/4323/](https://peerj.com/articles/4323/)
    - The overall most abundant, sOTU1 (39.5% average relative abundance) was classified to clade C, sOTU2 (21.7%) to clade A (A3) and sOTU3 (16.4%) and sOTU4 (10.2%) to clade B and type B1, respectively. sOTU5, classified as an H1 type, was the fifth most abundant sOTU (4.8%), and dominated in one of the 14 samples (Fig. 1A). sOTU6 was classified as a G2-1 type and was the sixth most abundant sOTU (3.3%). The one sOTU classified to clade D, type D1 (Symbiodinium trenchii, was rare and occupied the 24th position.
  - C1 *Cladocopium goreaui* [https://doi.org/10.3389/fmars.2018.00150](https://doi.org/10.3389/fmars.2018.00150)


----------------------------------------------------------------------------------------------
## Adapter trimming-clipping-quality filtering
### fastx_toolkit

Siderastrea *clippedtrimmed_nofilter.fastq and associated FASTQC files:
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/fastqc/

The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. 
FASTX Toolkit 0.0.13 by A. Gordon (gordon@cshl.edu)
[http://hannonlab.cshl.edu/fastx_toolkit/](http://hannonlab.cshl.edu/fastx_toolkit/)


This should be run from within the folder with all of your original .fastq files
Will quality trim multiple SINGLE-END fastq files

Things to customize for your particular platform:
- qualoffset = 33 Quality score offset
- -t option in qualitytrim (the lower threshold quality score for trimming):  -t 20 
- -l option in qualitytrim and adapterclip (the length threshold for throwing out short reads):  -l 20
- -q -p options in quality filter (the low quality, -q, and percentage -p options for filtration):  -q 20 -p 90

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq/originalfastqs/adapterlist.txt /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

```
nano TrimClipFilter.sh
```

```
#!/bin/bash -l

#SBATCH -o TrimClipFilter.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=TrimClipFilter

/cm/shared/courses/dbarshis/15AdvBioinf/scripts/Trimclipfilterstatsbatch_advbioinf.py adapterlist.txt *.fastq
```

```
sbatch TrimClipFilter.sh
```


----------------------------------------------------------------------------------------------
## QA-QC
### Examine quality metrics for sequencing reads
- Examining the base quality distribution, kmer frequencies and adapter contamination by position in the read is an important first step to understanding the underlying quality of your data. 
- For example, an increase in adapter frequency as one moves along a read is indicative of incomplete removal of adapter sequence during demultiplexing, a rather common occurrence. 
- In addition, the presence of over-represented sequences can be indicative of adapter contamination, rRNA reads, or perhaps other exongenous contamination.

Aims: 
- run trim clipped filtered stats Python script on stats output file 
- check adapter clipping stats for each file

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq/filteringstats
```
nano filterstats.sh
```

```
#!/bin/bash -l

#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=filter-stats

/cm/shared/courses/dbarshis/15AdvBioinf/scripts/Schafran_trimstatstable_advbioinf.py Sid_trimclipstats_final.txt Sid_trimclipstats_out.txt
```

```
sbatch filterstats.sh
```

### Filtering stats
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/filteringstats/

cp -avr /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq/filteringstats/ /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Original output TrimClipFilter.txt, renamed as Sid_trimclipstats_final.txt
```
cat Sid_trimclipstats_final.txt
```

Each individual sample has the following associated files:
- _clipped_trimmed_stats.txt
- _nucdist.png
- _qualbox.png


--------------------------------------------------------------------------------------------
## Test *Siderastrea siderea* Davies reference transcriptome

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

### *Siderastrea siderea* reference transcriptome
- File:  davies_Ssid.fasta
- from Assistant Professor Sarah Davies, Boston University [https://sites.bu.edu/davieslab/data-code/](https://sites.bu.edu/davieslab/data-code/)
- Assembled and annotated transcriptome for the adult scleractinian coral Siderastrea siderea with symbiont contamination removed. 
- Data were generated using Illumina HiSeq2000 2*100bp reads and assembled using Trinity. 
- Citaion: Davies SW, Marchetti A, Ries, JB and KD Castillo (2016). Thermal and pCO2 stress elicit divergent transcriptomic responses in a resilient coral. Frontiers in Marine Science. FMARS-03-00112.
[https://www.frontiersin.org/articles/10.3389/fmars.2016.00112/full](https://www.frontiersin.org/articles/10.3389/fmars.2016.00112/full)

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py davies_Ssid.fasta
```

```
The total number of sequences is 46704
The average sequence length is 2233
The total number of bases is 104309788
The minimum sequence length is 500
The maximum sequence length is 18031
The N50 is 2747
Median Length = 1499
contigs < 150bp = 0
contigs >= 500bp = 46704
contigs >= 1000bp = 38705
contigs >= 2000bp = 22065
```

**Rename fasta file based on # of contigs**
```
mv davies_Ssid.fasta 46704_davies_Ssid.fasta
```

### Check assembly with Trinity stats script

```
enable_lmod
module load container_env trinity
crun TrinityStats.pl 46704_davies_Ssid.fasta
```

Output:

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	17322
Total trinity transcripts:	46704
Percent GC: 42.99

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 5422
	Contig N20: 4370
	Contig N30: 3661
	Contig N40: 3161
	Contig N50: 2747

	Median contig length: 1915
	Average contig: 2233.42
	Total assembled bases: 104309788


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 5722
	Contig N20: 4529
	Contig N30: 3756
	Contig N40: 3211
	Contig N50: 2752

	Median contig length: 1790
	Average contig: 2159.32
	Total assembled bases: 37403676



### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### bowtie for Siderastrea mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/46704_davies_Ssid.fasta sid
```

## Test coral host *Siderastrea siderea* transcriptome
Mapping to reference
```
nano mapreads_sid_2.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_sid_2.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_sid_2

module load bowtie2/2.2.4

bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_349_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_349_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_366_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_366_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_376_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_376_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_331_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_331_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_338_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_338_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_348_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_348_N_Sid_nof_sid.sam -k 5\n
```

```
sbatch mapreads_sid_2.sh
```

```
cat bowtie2_sid_2.txt 
```

##### NOT FILTERED
```
10127834 reads; of these:
  10127834 (100.00%) were unpaired; of these:
    9554732 (94.34%) aligned 0 times
    401788 (3.97%) aligned exactly 1 time
    171314 (1.69%) aligned >1 times
5.66% overall alignment rate

19306925 reads; of these:
  19306925 (100.00%) were unpaired; of these:
    18389915 (95.25%) aligned 0 times
    692496 (3.59%) aligned exactly 1 time
    224514 (1.16%) aligned >1 times
4.75% overall alignment rate

12103752 reads; of these:
  12103752 (100.00%) were unpaired; of these:
    10850345 (89.64%) aligned 0 times
    381801 (3.15%) aligned exactly 1 time
    871606 (7.20%) aligned >1 times
10.36% overall alignment rate

8479518 reads; of these:
  8479518 (100.00%) were unpaired; of these:
    7986208 (94.18%) aligned 0 times
    246952 (2.91%) aligned exactly 1 time
    246358 (2.91%) aligned >1 times
5.82% overall alignment rate

4539491 reads; of these:
  4539491 (100.00%) were unpaired; of these:
    3813868 (84.02%) aligned 0 times
    121689 (2.68%) aligned exactly 1 time
    603934 (13.30%) aligned >1 times
15.98% overall alignment rate

23328804 reads; of these:
  23328804 (100.00%) were unpaired; of these:
    22049555 (94.52%) aligned 0 times
    900432 (3.86%) aligned exactly 1 time
    378817 (1.62%) aligned >1 times
5.48% overall alignment rate
```

## Outcome host mapping to Davies Siderastrea reference transcriptome
- tested 6 samples
- 5-16% overall alignment
- ~3-4% singly aligned reads
- This is low alignment


--------------------------------------------------------------------------------------------
## Test *Siderastrea siderea* Davies reference transcriptome - additional test mapping

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/test/

### Make transcriptome mappable
Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### bowtie for Siderastrea mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/46704_davies_Ssid.fasta sid
```

## Test coral host *Siderastrea siderea* reference transcriptome
Mapping to reference
```
nano mapreads_sid_test.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_sid_test.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_sid_test

module load bowtie2/2.2.4

bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_361_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_361_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_371_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_371_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_335_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_335_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_342_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_342_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_347_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_347_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_003_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_003_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_011_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_011_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_019_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_019_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_026_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_026_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_030_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_030_C_Sid_nof_sid.sam -k 5\n
```

```
sbatch mapreads_sid_test.sh
```

```
cat bowtie2_sid_test.txt 
```

##### NOT FILTERED
```
12959491 reads; of these:
  12959491 (100.00%) were unpaired; of these:
    12303149 (94.94%) aligned 0 times
    396652 (3.06%) aligned exactly 1 time
    259690 (2.00%) aligned >1 times
5.06% overall alignment rate

18367528 reads; of these:
  18367528 (100.00%) were unpaired; of these:
    17321445 (94.30%) aligned 0 times
    792353 (4.31%) aligned exactly 1 time
    253730 (1.38%) aligned >1 times
5.70% overall alignment rate

8523378 reads; of these:
  8523378 (100.00%) were unpaired; of these:
    8091660 (94.93%) aligned 0 times
    267455 (3.14%) aligned exactly 1 time
    164263 (1.93%) aligned >1 times
5.07% overall alignment rate

15309393 reads; of these:
  15309393 (100.00%) were unpaired; of these:
    14504424 (94.74%) aligned 0 times
    516896 (3.38%) aligned exactly 1 time
    288073 (1.88%) aligned >1 times
5.26% overall alignment rate

25535025 reads; of these:
  25535025 (100.00%) were unpaired; of these:
    24360145 (95.40%) aligned 0 times
    872482 (3.42%) aligned exactly 1 time
    302398 (1.18%) aligned >1 times
4.60% overall alignment rate

13100851 reads; of these:
  13100851 (100.00%) were unpaired; of these:
    12505130 (95.45%) aligned 0 times
    400055 (3.05%) aligned exactly 1 time
    195666 (1.49%) aligned >1 times
4.55% overall alignment rate

16886305 reads; of these:
  16886305 (100.00%) were unpaired; of these:
    16062593 (95.12%) aligned 0 times
    599966 (3.55%) aligned exactly 1 time
    223746 (1.33%) aligned >1 times
4.88% overall alignment rate

23501378 reads; of these:
  23501378 (100.00%) were unpaired; of these:
    22524930 (95.85%) aligned 0 times
    736595 (3.13%) aligned exactly 1 time
    239853 (1.02%) aligned >1 times
4.15% overall alignment rate

12065549 reads; of these:
  12065549 (100.00%) were unpaired; of these:
    10807842 (89.58%) aligned 0 times
    389263 (3.23%) aligned exactly 1 time
    868444 (7.20%) aligned >1 times
10.42% overall alignment rate

21883993 reads; of these:
  21883993 (100.00%) were unpaired; of these:
    21011965 (96.02%) aligned 0 times
    642534 (2.94%) aligned exactly 1 time
    229494 (1.05%) aligned >1 times
3.98% overall alignment rate
```

## Outcome host mapping to Davies Siderastrea reference
- tested additional 10 samples
- 4-10% overall alignment
- 3-4% singly aligned

***Very low alignment, need to investigate further***


--------------------------------------------------------------------------------------------
## Test *Siderastrea siderea* Davies reference transcriptome - additional test mapping Siderastrea samples with highest retained reads
- n=9 samples with >=85% retained reads when filtered or not filtered

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/test/

### Make transcriptome mappable
Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### bowtie for Siderastrea mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/46704_davies_Ssid.fasta sid
```

## Test coral host *Siderastrea siderea* reference transcriptome
Mapping to reference
```
nano mapreads_sid_test_2.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_sid_test_2.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_sid_test_2

module load bowtie2/2.2.4

bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_012_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_012_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_361_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_361_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_364_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_364_N_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_334_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_334_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_022_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_022_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_024_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_024_C_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_341_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_341_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_004_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_004_La_Sid_nof_sid.sam -k 5\n
bowtie2 --local -x sid -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_013_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_013_La_Sid_nof_sid.sam -k 5\n
```

```
sbatch mapreads_sid_test_2.sh
```

```
cat bowtie2_sid_test_2.txt 
```

##### NOT FILTERED
```
17696748 reads; of these:
  17696748 (100.00%) were unpaired; of these:
    15569908 (87.98%) aligned 0 times
    1625553 (9.19%) aligned exactly 1 time
    501287 (2.83%) aligned >1 times
12.02% overall alignment rate

12959491 reads; of these:
  12959491 (100.00%) were unpaired; of these:
    12303149 (94.94%) aligned 0 times
    396652 (3.06%) aligned exactly 1 time
    259690 (2.00%) aligned >1 times
5.06% overall alignment rate

22723159 reads; of these:
  22723159 (100.00%) were unpaired; of these:
    21661238 (95.33%) aligned 0 times
    792810 (3.49%) aligned exactly 1 time
    269111 (1.18%) aligned >1 times
4.67% overall alignment rate

18032416 reads; of these:
  18032416 (100.00%) were unpaired; of these:
    16974956 (94.14%) aligned 0 times
    680487 (3.77%) aligned exactly 1 time
    376973 (2.09%) aligned >1 times
5.86% overall alignment rate

14496181 reads; of these:
  14496181 (100.00%) were unpaired; of these:
    13377832 (92.29%) aligned 0 times
    904733 (6.24%) aligned exactly 1 time
    213616 (1.47%) aligned >1 times
7.71% overall alignment rate

16287244 reads; of these:
  16287244 (100.00%) were unpaired; of these:
    15138875 (92.95%) aligned 0 times
    845760 (5.19%) aligned exactly 1 time
    302609 (1.86%) aligned >1 times
7.05% overall alignment rate

17479996 reads; of these:
  17479996 (100.00%) were unpaired; of these:
    16683471 (95.44%) aligned 0 times
    596472 (3.41%) aligned exactly 1 time
    200053 (1.14%) aligned >1 times
4.56% overall alignment rate

20666252 reads; of these:
  20666252 (100.00%) were unpaired; of these:
    19681506 (95.24%) aligned 0 times
    722461 (3.50%) aligned exactly 1 time
    262285 (1.27%) aligned >1 times
4.76% overall alignment rate

20281418 reads; of these:
  20281418 (100.00%) were unpaired; of these:
    18627958 (91.85%) aligned 0 times
    1102462 (5.44%) aligned exactly 1 time
    550998 (2.72%) aligned >1 times
8.15% overall alignment rate
```

## Outcome host mapping to Davies Siderastrea transcriptome reference
- tested additional 9 samples
- 4-12% overall alignment
- 3-9% singly aligned

***Very low alignment, need to investigate further***


--------------------------------------------------------------------------------------------
## Test *Siderastrea radians* host reference assembly
- from Dr. Viridiana Avila-Magaña & Professor Mónica Medina, Pennsylvania State University 
[http://medinalab.org/new/](http://medinalab.org/new/)
- contains the host transcripts retrieved by a blast against coral, Aiptasia and Nematostella genomes
- File: Sid_Host.fna

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Sid_Host.fna
```

```
The total number of sequences is 685205
The average sequence length is 585
The total number of bases is 401076233
The minimum sequence length is 201
The maximum sequence length is 2459
The N50 is 779
Median Length = 1025
contigs < 150bp = 0
contigs >= 500bp = 260979
contigs >= 1000bp = 105374
contigs >= 2000bp = 17025
```

Rename fasta file based on # of contigs
```
mv Sid_Host.fna 685205_Sid_radians_Avila-Medina.fna
```

### Check assembly
```
enable_lmod
module load container_env trinity
crun TrinityStats.pl 685205_Sid_radians_Avila-Medina.fna
```

Output:

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	475058
Total trinity transcripts:	685205
Percent GC: 40.37

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 1973
	Contig N20: 1589
	Contig N30: 1267
	Contig N40: 1002
	Contig N50: 779

	Median contig length: 392
	Average contig: 585.34
	Total assembled bases: 401076233


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 1951
	Contig N20: 1423
	Contig N30: 1013
	Contig N40: 752
	Contig N50: 574

	Median contig length: 335
	Average contig: 494.07
	Total assembled bases: 234713656


### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/685205_Sid_radians_Avila-Medina.fna Sid_radians
```

## Test mapping to coral host *Siderastrea radians* reference transcriptome
- test mapping our Siderastrea samples with the highest retained reads
- n=14 samples with >=82% retained reads when filtered or not filtered

Mapping to reference
```
nano mapreads_Sid_radians.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Sid_radians.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Sid_radians

module load bowtie2/2.2.4

bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_025_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_025_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_372_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_372_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_019_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_019_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_337_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_337_C_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_368_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_368_C_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_012_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_012_C_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_361_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_361_N_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PC_364_N_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PC_364_N_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_334_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_334_C_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_022_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_022_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_024_C_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_024_C_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/PO_341_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S PO_341_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_004_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_004_La_Sid_nof_Sid_radians.sam -k 5\n
bowtie2 --local -x Sid_radians -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/R_013_La_Sid_yr2_R1_clippedtrimmed_nofilter.fastq -S R_013_La_Sid_nof_Sid_radians.sam -k 5\n
```

```
sbatch mapreads_Sid_radians.sh
```

```
cat bowtie2_Sid_radians.txt 
```

```
23501378 reads; of these:
  23501378 (100.00%) were unpaired; of these:
    2106172 (8.96%) aligned 0 times
    1303159 (5.55%) aligned exactly 1 time
    20092047 (85.49%) aligned >1 times
91.04% overall alignment rate

12959491 reads; of these:
  12959491 (100.00%) were unpaired; of these:
    2284328 (17.63%) aligned 0 times
    616226 (4.76%) aligned exactly 1 time
    10058937 (77.62%) aligned >1 times
82.37% overall alignment rate

17696748 reads; of these:
  17696748 (100.00%) were unpaired; of these:
    4824605 (27.26%) aligned 0 times
    1183624 (6.69%) aligned exactly 1 time
    11688519 (66.05%) aligned >1 times
72.74% overall alignment rate

22723159 reads; of these:
  22723159 (100.00%) were unpaired; of these:
    1989520 (8.76%) aligned 0 times
    1109852 (4.88%) aligned exactly 1 time
    19623787 (86.36%) aligned >1 times
91.24% overall alignment rate

18032416 reads; of these:
  18032416 (100.00%) were unpaired; of these:
    3415945 (18.94%) aligned 0 times
    995994 (5.52%) aligned exactly 1 time
    13620477 (75.53%) aligned >1 times
81.06% overall alignment rate

14496181 reads; of these:
  14496181 (100.00%) were unpaired; of these:
    2203207 (15.20%) aligned 0 times
    808110 (5.57%) aligned exactly 1 time
    11484864 (79.23%) aligned >1 times
84.80% overall alignment rate

16287244 reads; of these:
  16287244 (100.00%) were unpaired; of these:
    3451586 (21.19%) aligned 0 times
    865214 (5.31%) aligned exactly 1 time
    11970444 (73.50%) aligned >1 times
78.81% overall alignment rate

17479996 reads; of these:
  17479996 (100.00%) were unpaired; of these:
    1891835 (10.82%) aligned 0 times
    935835 (5.35%) aligned exactly 1 time
    14652326 (83.82%) aligned >1 times
89.18% overall alignment rate

20666252 reads; of these:
  20666252 (100.00%) were unpaired; of these:
    2000551 (9.68%) aligned 0 times
    1180083 (5.71%) aligned exactly 1 time
    17485618 (84.61%) aligned >1 times
90.32% overall alignment rate

20281418 reads; of these:
  20281418 (100.00%) were unpaired; of these:
    3543869 (17.47%) aligned 0 times
    1359030 (6.70%) aligned exactly 1 time
    15378519 (75.83%) aligned >1 times
82.53% overall alignment rate

19967419 reads; of these:
  19967419 (100.00%) were unpaired; of these:
    1792457 (8.98%) aligned 0 times
    1084079 (5.43%) aligned exactly 1 time
    17090883 (85.59%) aligned >1 times
91.02% overall alignment rate

21243157 reads; of these:
  21243157 (100.00%) were unpaired; of these:
    2464607 (11.60%) aligned 0 times
    994128 (4.68%) aligned exactly 1 time
    17784422 (83.72%) aligned >1 times
88.40% overall alignment rate

16428742 reads; of these:
  16428742 (100.00%) were unpaired; of these:
    5297079 (32.24%) aligned 0 times
    1401396 (8.53%) aligned exactly 1 time
    9730267 (59.23%) aligned >1 times
67.76% overall alignment rate

22832214 reads; of these:
  22832214 (100.00%) were unpaired; of these:
    2890385 (12.66%) aligned 0 times
    904271 (3.96%) aligned exactly 1 time
    19037558 (83.38%) aligned >1 times
87.34% overall alignment rate
```

Outcome:
- 67-91% overall alignment
- 4-9% singly aligned


--------------------------------------------------------------------------------------------
## Decision about Host reference

Still very low alignment to our specific population.
**I did a de novo transcriptome assembly of our specific popluation. This will be used to make a hybrid reference assembly below.**

--------------------------------------------------------------------------------------------
## Identify dominant symbiont type

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/symbiont_typing/

- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- map samples against ITS2 database
- Arif et al. 2014 [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
- Corrigendum, with updated ITS2 database file, published 17 July 2019 [https://onlinelibrary.wiley.com/doi/10.1111/mec.14956](https://onlinelibrary.wiley.com/doi/10.1111/mec.14956)
- File:  mec14956-sup-0001-files1_corrigendum.fasta

***use revised Arif ITS2 database***
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py mec14956-sup-0001-files1_corrigendum.fasta
```

```
The total number of sequences is 400
The average sequence length is 376
The total number of bases is 150400
The minimum sequence length is 376
The maximum sequence length is 376
The N50 is 376
Median Length = 376
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------
## Map all samples to Arif ITS2 (all clades) database

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

### Make file mappable
Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

Bowtie files for Arif ITS2 clade mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mec14956-sup-0001-files1_corrigendum.fasta Arif_ITS2_corrigendum
```

#### Mapping to Arif_ITS2_corrigendum reference
```
nano mapreads_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Arif_ITS2_corrigendum

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Arif_ITS2_corrigendum -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Arif_ITS2_corrigendum.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2_corrigendum.sh
```

Output:  bowtie2_Arif_ITS2_corrigendum.txt
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/symbiont_typing/


--------------------------------------------------------------------------------------------
## Count expression - all reads mapped to Arif_ITS2 symbiont database reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

```
nano countexpression_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_Arif_ITS2_corrigendum_Sid

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_Arif_ITS2_corrigendum.sam
```

```
sbatch countexpression_Arif_ITS2_corrigendum.sh
```

*Output match_counts.txt file renamed as Siderastrea_symbiont-references_mapping_match-counts_final.xlsx in local folder on computer.*


--------------------------------------------------------------------------------------------
## Merge all *nof_Arif_ITS2_counts.txt files into one big table 

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/symbiont_typing/

- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 


```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=append_filename

for f in *_Arif_ITS2_corrigendum_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

#### Concatenate files
```
cat *_nof_Arif_ITS2_corrigendum_counts.txt > merged_Sid_ITS2_counts_Arif-corrigendum.txt
```

*Copied file merged_Sid_ITS2_counts_Arif-corrigendum.txt to local folder on computer*
File:  merged_Sid_ITS2_counts_Arif-corrigendum.xlsx

## Outcome symbiont clade mapping
- 90% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 72% mapped to C
  - 15% mapped to B
- Need to test the few different C and B references
- Then pick the two with the best mapping rate, merge them 
- Then re-map to see how much we lose to multiply mapping contigs

File:  Siderastrea symbiont refs mapping.xlsx

--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome - C1 *Cladocopium goreaui* Davies reference
- Davies et al. (2018). Symbiodinium functional diversity in the coral Siderastrea siderea is influenced by thermal stress and reef environment, but not ocean acidification. Frontiers in Marine Science. FMARS-05-00150.
- Assembled and annotated transcriptome for the symbiotic dinoflagellate algae Cladocopium hosted by *Siderastrea siderea* with all host contamination removed. 
- Data were generated using Illumina HiSeq2000 2*100bp reads and assembled using Trinity.
- File:  davies_cladeC_feb.fasta
[https://sites.bu.edu/davieslab/data-code/](https://sites.bu.edu/davieslab/data-code/)
[https://doi.org/10.3389/fmars.2018.00150](https://doi.org/10.3389/fmars.2018.00150)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py davies_cladeC_feb.fasta
```

```
The total number of sequences is 65838
The average sequence length is 1482
The total number of bases is 97581498
The minimum sequence length is 500
The maximum sequence length is 18168
The N50 is 1746
Median Length = 1297
contigs < 150bp = 0
contigs >= 500bp = 65838
contigs >= 1000bp = 40840
contigs >= 2000bp = 13246
```

### Rename fasta file based on # of contigs
```
mv davies_cladeC_feb.fasta 65838_davies_cladeC_feb.fasta
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb.fasta davies_cladeC
```

## Test Davies C1 *Cladocopium goreaui* reference
Mapping to reference
```
nano mapreads_davies_cladeC.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_davies_cladeC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_davies_cladeC

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x davies_cladeC -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_davies_cladeC.sam -k 5\n; done
```

```
sbatch mapreads_davies_cladeC.sh
```

```
cat bowtie2_davies_cladeC.txt 
```

Outcome:
- 2-10% overall alignment (avg. 4.4%)
- 1.5-6.5% singly aligned (avg. 3.2%)

## Count expression 
```
nano countexpression_davies_cladeC.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_davies_cladeC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_davies_cladeC

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_davies_cladeC.sam
```

```
sbatch countexpression_davies_cladeC.sh
```

--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome -  *Cladocopium goreaui* reference (ITS2 type C1)
- Chen et al. (2020). Evidence That Inconsistent Gene Prediction Can Mislead Analysis of Dinoflagellate Genomes. Journal of Phycology.
- predicted protein‐coding genes from published draft Symbiodiniaceae genome of *Cladocopium goreaui* from Liu et al. (2018)
    - Liu et al. (2018). Symbiodinium genomes reveal adaptive evolution of functions related to coral‐dinoflagellate symbiosis. Commun. Biol. 1:95.
        - Symbiodinium goreaui (Clade C, type C1; AIMS-aten-C1-MI-cfu-B2, now AIMS culture collection SCF055-01) is a single-cell monoclonal culture first isolated from the coral Acropora tenuis at Magnetic Island (Queensland, Australia) at 3 m depth; this culture is maintained at the Australian Institute of Marine Science, Townsville, Australia.
        - Illumina HiSeq 2500 platform, 2 × 150 bp reads
    - generated 116.0 Gb (614.6 million reads)
- file:  Cladocopium_goreaui.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
- Chen et al. (2019). Revised genome sequences and annotations of six Symbiodiniaceae taxa.
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Cladocopium_goreaui.CDS.fna
```

```
The total number of sequences is 39006
The average sequence length is 1625
The total number of bases is 63419290
The minimum sequence length is 111
The maximum sequence length is 40317
The N50 is 2388
Median Length = 852
contigs < 150bp = 5
contigs >= 500bp = 31703
contigs >= 1000bp = 21315
contigs >= 2000bp = 9624
```

#### Rename fasta file based on # of contigs
```
mv Cladocopium_goreaui.CDS.fna 39006_Cladocopium_goreaui.CDS.fna
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/39006_Cladocopium_goreaui.CDS.fna Cladocopium_goreaui
```

## Test Cladocopium_goreaui reference
Mapping to reference
```
nano mapreads_Cladocopium_goreaui.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium_goreaui.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Cladocopium_goreaui

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq; do
        bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
                     --rg SM:${i%_clippedtrimmed_nofilter.fastq} \
                     --local -x Cladocopium_goreaui -U $i \
                        > ${i%_clippedtrimmed_nofilter.fastq}_nof_Cladocopium_goreaui.sam -k 5\n;
done
```

```
sbatch mapreads_Cladocopium_goreaui.sh
```

```
cat bowtie2_Cladocopium_goreaui.txt 
```

Outcome
- 0-4.8% overall alignment
- 0-4% singly aligned

## Count expression 
```
nano countexpression_Cladocopium_goreaui.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Cladocopium_goreaui.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Cladocopium_goreaui

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_Cladocopium_goreaui.sam
```

```
sbatch countexpression_Cladocopium_goreaui.sh
```


--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome -  *Cladocopium sp. C92* reference
- Chen et al. (2020). Evidence That Inconsistent Gene Prediction Can Mislead Analysis of Dinoflagellate Genomes. Journal of Phycology.
- predicted protein‐coding genes from published draft Symbiodiniaceae genome of *Cladocopium C92* from Shoguchi et al. (2018)
    - Shoguchi et al. (2018). Two divergent Symbiodinium genomes reveal conservation of a gene cluster for sunscreen biosynthesis and recently lost genes. BMC Genomics 19:458.
        - Dinoflagellate Symbiodinium spp. clade C (SymC) were cultured to produce genomic DNA and mRNA for sequencing. SymC are harbored by the cardiid clams Fragum sp. obtained in Okinawa, Japan. In regard to host habitats, Fragum is infaunal [72]. In the 1980s, isolations of Symbiodinium cells were performed by Prof. Terufumi Yamasu at the University of the Ryukyus using sterilized seawater and micropipettes. The cultured Symbiodinium have been maintained since then in the laboratory of Prof. Michio Hidaka, at the University of the Ryukyus. SymC were designated as strain “Y103”. By manually isolating single cells under a microscope using a glass micropipette, each isoclonal line was established at the Marine Genomics Unit of Okinawa Institute of Science and Technology Graduate University in 2009. Repetitive subculture in 250-mL flasks has continued for 8 years.
        - Illumina Genome Analyzer IIx (GAIIx) and Hiseq Paired-end reads
- file:  Cladocopium_sp_C92.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
- Chen et al. (2019). Revised genome sequences and annotations of six Symbiodiniaceae taxa.
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Cladocopium_sp_C92.CDS.fna
```

```
The total number of sequences is 33421
The average sequence length is 2201
The total number of bases is 73570391
The minimum sequence length is 93
The maximum sequence length is 28947
The N50 is 3198
Median Length = 2229
contigs < 150bp = 1
contigs >= 500bp = 30177
contigs >= 1000bp = 23167
contigs >= 2000bp = 12664
```

#### Rename fasta file based on # of contigs
```
mv Cladocopium_sp_C92.CDS.fna 33421_Cladocopium_sp_C92.CDS.fna
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/33421_Cladocopium_sp_C92.CDS.fna Cladocopium_sp_C92
```

## Test Cladocopium_sp_C92 reference
Mapping to reference
```
nano mapreads_Cladocopium_sp_C92.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium_sp_C92.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Cladocopium_sp_C92

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Cladocopium_sp_C92 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Cladocopium_sp_C92.sam -k 5\n; done
```

```
sbatch mapreads_Cladocopium_sp_C92.sh
```

```
cat bowtie2_Cladocopium_sp_C92.txt 
```

Outcome
- 0-3.7% overall alignment
- ~0-2.5% singly aligned

## Count expression 
```
nano countexpression_Cladocopium_sp_C92.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Cladocopium_sp_C92.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Cladocopium_sp_C92

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_Cladocopium_sp_C92.sam
```

```
sbatch countexpression_Cladocopium_sp_C92.sh
```


--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome -  *Breviolum minutum* (ITS2 type B1) reference
- Chen et al. (2020). Evidence That Inconsistent Gene Prediction Can Mislead Analysis of Dinoflagellate Genomes. Journal of Phycology.
- predicted protein‐coding genes from published draft Symbiodiniaceae genome of *Breviolum minutum* from Shoguchi et al. (2013)
    - Shoguchi et al. (2013). Draft Assembly of the Symbiodinium minutum Nuclear Genome Reveals Dinoflagellate Gene Structure. Current Biology.
        - DNA obtained from a single clonal culture of the Symbiodinium minutum was sequenced with Roche 454 GS-FLX and the Illumina Genome Analyzer IIx (GAIIx). 
- file:  Breviolum_minutum.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
- Chen et al. (2019). Revised genome sequences and annotations of six Symbiodiniaceae taxa.
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Breviolum_minutum.CDS.fna
```

```
The total number of sequences is 32803
The average sequence length is 1929
The total number of bases is 63302208
The minimum sequence length is 102
The maximum sequence length is 32739
The N50 is 2772
Median Length = 903
contigs < 150bp = 2
contigs >= 500bp = 28536
contigs >= 1000bp = 20813
contigs >= 2000bp = 10169
```

#### Rename fasta file based on # of contigs
```
mv Breviolum_minutum.CDS.fna 32803_Breviolum_minutum.CDS.fna
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/32803_Breviolum_minutum.CDS.fna Breviolum_minutum
```

## Test Breviolum_minutum reference
Mapping to reference
```
nano mapreads_Breviolum_minutum.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Breviolum_minutum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Breviolum_minutum

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Breviolum_minutum -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Breviolum_minutum.sam -k 5\n; done
```

```
sbatch mapreads_Breviolum_minutum.sh
```

```
cat bowtie2_Breviolum_minutum.txt 
```

Outcome:
- 0.8% overall alignment
- 0.6% singly aligned

## Count expression 
```
nano countexpression_Breviolum_minutum.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Breviolum_minutum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Breviolum_minutum

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_Breviolum_minutum.sam
```

```
sbatch countexpression_Breviolum_minutum.sh
```


--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome -  *Breviolum psygmophilum* (ITS2 type B2) reference
- Davies et al. (2018). Symbiodinium Functional Diversity in the Coral Siderastrea siderea Is Influenced by Thermal Stress and Reef Environment, but Not Ocean Acidification. Frontiers.
- file:  B_psygmophilum_transcriptome.fasta
[https://www.frontiersin.org/articles/10.3389/fmars.2018.00150/full](https://www.frontiersin.org/articles/10.3389/fmars.2018.00150/full)
[http://sites.bu.edu/davieslab/data-code/](http://sites.bu.edu/davieslab/data-code/)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py B_psygmophilum_transcriptome.fasta
```

```
The total number of sequences is 31970
The average sequence length is 1290
The total number of bases is 41260850
The minimum sequence length is 500
The maximum sequence length is 8436
The N50 is 1458
Median Length = 823
contigs < 150bp = 0
contigs >= 500bp = 31970
contigs >= 1000bp = 18541
contigs >= 2000bp = 4186
```

#### Rename fasta file based on # of contigs
```
mv B_psygmophilum_transcriptome.fasta 31970_Breviolum_psygmophilum_transcriptome.fasta
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/31970_Breviolum_psygmophilum_transcriptome.fasta Breviolum_psygmophilum
```

## Test Breviolum psygmophilum reference
Mapping to reference
```
nano mapreads_Breviolum_psygmophilum.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Breviolum_psygmophilum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Breviolum_psygmophilum

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Breviolum_psygmophilum -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Breviolum_psygmophilum.sam -k 5\n; done
```

```
sbatch mapreads_Breviolum_psygmophilum.sh
```

```
cat bowtie2_Breviolum_psygmophilum.txt 
```

Outcome:
- 0.7% overall alignment
- 0.5% singly aligned

## Count expression 
```
nano countexpression_Breviolum_psygmophilum.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Breviolum_psygmophilum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Breviolum_psygmophilum

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_Breviolum_psygmophilum.sam
```

```
sbatch countexpression_Breviolum_psygmophilum.sh
```


--------------------------------------------------------------------------------------------
## Test symbiont reference transcriptome -  *Breviolum B5* reference Medina-Avila
- from Dr. Viridiana Avila-Magaña & Professor Mónica Medina, Pennsylvania State University 
[http://medinalab.org/new/](http://medinalab.org/new/)
- isolated from *Siderastrea radians* that typically has B5 (while *S. siderea* typically hosts C1) in the Mexican Caribbean
- contains the symbiont transcripts retrieved by a blast against Symbiodiniaceae genomes
    - blasted the metatranscriptomes to the published Symbiodiniaceae genomes at that time (S. microadriaticum, B. minutum, and F. kawagutii). 
- They did Symbiodiniaceae ITS2 amplicon sequencing confirming that the S. radians symbiont identity was Breviolum (100%), and did not detect other background populations in these samples.
- File:  Symbio_Sider.fna

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/


## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Symbio_Sider.fna
```

```
The total number of sequences is 173968
The average sequence length is 627
The total number of bases is 109245170
The minimum sequence length is 201
The maximum sequence length is 2459
The N50 is 1040
Median Length = 258
contigs < 150bp = 0
contigs >= 500bp = 66248
contigs >= 1000bp = 36581
contigs >= 2000bp = 5807
```

#### Rename fasta file based on # of contigs
```
mv Symbio_Sider.fna 173968_Breviolum_B5_Sid_radians_Avila-Medina.fna
```

### Make transcriptome mappable
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/

Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

#### Bowtie for symbiont mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/173968_Breviolum_B5_Sid_radians_Avila-Medina.fna Breviolum_B5
```

## Test Breviolum B5 reference
Mapping to reference
```
nano mapreads_Breviolum_B5.sh
```

**Forgot to edit the associated name of the ___.sam file name in the script below, and therefore the Breviolum_B5 counts and sam files are actually labeled as _nof_Breviolum_psygmophilum_counts.txt and _nof_Breviolum_psygmophilum.sam**

```
#!/bin/bash -l

#SBATCH -o bowtie2_Breviolum_B5.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Breviolum_B5

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Breviolum_B5 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Breviolum_psygmophilum.sam -k 5\n; done
```

```
sbatch mapreads_Breviolum_B5.sh
```

```
cat bowtie2_Breviolum_B5.txt 
```

Outcome
- 13-49% overall alignment
- 4-17% singly aligned

## Count expression 
```
nano countexpression_Breviolum_B5.sh
```

*Should have been files *_nof_Breviolum_B5.sam but I made mistake in above script, and it overwrote the actual _nof_Breviolum_psygmophilum.sam files*
```
#!/bin/bash -l

#SBATCH -o countexpression_Breviolum_B5.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_Breviolum_B5

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_nof_Breviolum_psygmophilum.sam
```

```
sbatch countexpression_Breviolum_B5.sh
```


--------------------------------------------------------------------------------------------
## Decision about Symbiont reference

**Will use Davies clade C and Avila-Medina clade B symbiont transcriptome references**

--------------------------------------------------------------------------------------------
## Create hybrid reference for mapping

***This was the first mapping attempt, after which we realized that Davies and Avila-Medina references were not filtered the same way as our de novo assemblies. Also, the majority of our ~3000 de novo symbiont sequences overlapped with the clade C reference. So we are not using our symbiont de novo assembly.***
*See below for final hybrid ref used*

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/

Merge renamed Good Coral and Good Symbiont assembly with the relevant predicted transcripts from the Symbiont references to create a hybrid reference for mapping.

## Breviolum transcriptome reference

### Add suffix to fasta reference sequence names
addsuffixtofastaseqnames.py

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

```
nano addsuffixtofastaseqnames_Breviolum-B5.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_Breviolum-B5.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Symbiont /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/173968_Breviolum_B5_Sid_radians_Avila-Medina.fna
```

```
sbatch addsuffixtofastaseqnames_Breviolum-B5.sh
```

Output:   173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta

```
cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/173968_Breviolum_B5_Sid_radians_Avila-Medina_suffixed.fasta
```

## Cladocopium transcriptome reference

*Potential problem - each line of (same) sequence is on a new line?*
head -20 65838_davies_cladeC_feb.fasta
```
>comp3_c0_seq1 len=529 path=[1253:0-528]
TGGACGTGCGGCACCCTCTCCAATCGACACGAAGATGTGCGCAGCTGCACCCTGCTGGCC
CGGCTGTGCAACACCCGGTAGGGCACACACCTCGCGTGCTGCAGAACCACCGGCACCGGT
ACTCTCTGGCGGCAGCAGGGCTTCCTTGAGACGTTCCTCGGCCTCGAGCAGGTCCACGAA
CTGTCCGCCCACCTTGATGAGCATGTCGCTGCGCCCTCGAAACTCCAGCGACCGCGCCCC
AGAGCTGAGGGTCACCAGGTCGGAGGTCTTGAAGAAGGGCGCGCCACCATCCGTTGGGGG
ATGAACGGACGAAGAGGCCCCACGGTACCCGGTGGTCACCATGGGCCCACGCAAGCACAG
CTCTCCAGTACCGTCCTTGCAAAGCTCTCCATCTCGCAGCACCGCCACCTCGGCGCCCTC
CACCGCCGTGAAGACGCTGCGCCCGCTCCGCGCAGAGCGGCCCTCGCTGTAGAGGCTCAA
CCAGTACTCGGTACTGATCAGCAATTCCACCACTCTCCTCCTCTCCGAG
>comp12_c0_seq1 len=928 path=[1062:0-927]
AGTCAAACACCACCACTGGAAGGTTCCCACAAAGGTCATTCCTGGATTGGTCGGGTAGTC
CGGATCCTGACTCTCCCGAGCATCCAGGACCAGAATGTCACCACGCGTCATGAGACGATC
CTCCGTGGCCAAGTAGGCGAAGACTGGAGATCGCCGCACCTTGACAGTGACCGAGGCATA
TGCGGCATTGTTGGGGTTGTTGATGGCAAAGGTCTCCACCGTGAAGTTGTAAATGTTCAG
ATCAGCTCCGGAGCCAATGGGTTCCAGAACGAAGGGCAGAATCACCAAAGACCTGGTGGT
GTAAACAATGTCTGGCATGGACGACAAATCCAGTTGGCCGGTGGTCTCCAACCAACGATA
GCCCAACTGTCTGGACACTGTTTGACACTCCGAAGGAACACCATTGGCCAACAGCGACAC
TTGCTCTGTGCGATTTACCAAGAGCTCTTTGGGTCCCAAGATACTGACCATCGGTGCCGG
GAAGTCCAACTTGGTGATATCAACCAACTCGCTCTGAGAGAGGCCCCATCGACTGGTGAC
```

grep -c ">" 65838_davies_cladeC_feb.fasta
65838

wc -l 65838_davies_cladeC_feb.fasta
1724567

Removed len= and path= from this file ('cleaned').

*Fixed-cleaned file:  65838_davies_cladeC_feb_fixed-cleaned.txt*

Add suffix _Symbiont
```
nano addsuffixtofastaseqnames_davies_cladeC.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_davies_cladeC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Symbiont /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-cleaned.txt
```

```
sbatch addsuffixtofastaseqnames_davies_cladeC.sh
```

Output:   65838_davies_cladeC_feb_fixed-clean_suffixed.fasta

```
cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/
```

#### Copy de novo assemblies
```
cp /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/3081_Sid_GoodSymb_Final_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/

cp /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/19272_Sid_GoodCoral_Final_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/
```

### Symbiont
*De novo assembly*
3081_Sid_GoodSymb_Final_renamed.fasta

173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta

65838_davies_cladeC_feb.fasta

### Host
*De novo assembly*
19272_Sid_GoodCoral_Final_renamed.fasta


--------------------------------------------------------------------------------------------
## Concatenate host reference and symbiont references to make mapping file

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/19272_Sid_GoodCoral_Final_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/3081_Sid_GoodSymb_Final_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/173968_Breviolum_B5_Sid_radians_Avila-Medina_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta > Sid_hybridreference.fasta
```

## Check assembly details

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Sid_hybridreference.fasta
```

```
The total number of sequences is 262159
The average sequence length is 881
The total number of bases is 231114187
The minimum sequence length is 201
The maximum sequence length is 32630
The N50 is 1342
Median Length = 244
contigs < 150bp = 0
contigs >= 500bp = 154372
contigs >= 1000bp = 84771
contigs >= 2000bp = 20963
```

--------------------------------------------------------------------------------------------
## Map all the samples to concatenated (host plus symbiont) hybrid reference transcriptome

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/Sid_hybridreference.fasta hybridref
```

#### Mapping to hybridreference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/

```
nano mapreads_hybridref.sh
```

For files *_Sid_yr2_R1_clippedtrimmed_nofilter.fastq
```
#!/bin/bash -l

#SBATCH -o bowtie2_hybridref.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_hybridref

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x hybridref -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_hybridref.sam -k 5\n; done
```

```
sbatch mapreads_hybridref.sh
```

```
cat bowtie2_hybridref.txt
```

```
10127834 reads; of these:
  10127834 (100.00%) were unpaired; of these:
    4602332 (45.44%) aligned 0 times
    2068383 (20.42%) aligned exactly 1 time
    3457119 (34.13%) aligned >1 times
54.56% overall alignment rate

4535297 reads; of these:
  4535297 (100.00%) were unpaired; of these:
    1873847 (41.32%) aligned 0 times
    930687 (20.52%) aligned exactly 1 time
    1730763 (38.16%) aligned >1 times
58.68% overall alignment rate

12959491 reads; of these:
  12959491 (100.00%) were unpaired; of these:
    6635412 (51.20%) aligned 0 times
    2552105 (19.69%) aligned exactly 1 time
    3771974 (29.11%) aligned >1 times
48.80% overall alignment rate

14417329 reads; of these:
  14417329 (100.00%) were unpaired; of these:
    6574097 (45.60%) aligned 0 times
    3077839 (21.35%) aligned exactly 1 time
    4765393 (33.05%) aligned >1 times
54.40% overall alignment rate

15619667 reads; of these:
  15619667 (100.00%) were unpaired; of these:
    7026683 (44.99%) aligned 0 times
    3253106 (20.83%) aligned exactly 1 time
    5339878 (34.19%) aligned >1 times
55.01% overall alignment rate

20103321 reads; of these:
  20103321 (100.00%) were unpaired; of these:
    14994975 (74.59%) aligned 0 times
    2227243 (11.08%) aligned exactly 1 time
    2881103 (14.33%) aligned >1 times
25.41% overall alignment rate

18367528 reads; of these:
  18367528 (100.00%) were unpaired; of these:
    8312671 (45.26%) aligned 0 times
    3631483 (19.77%) aligned exactly 1 time
    6423374 (34.97%) aligned >1 times
54.74% overall alignment rate

12449859 reads; of these:
  12449859 (100.00%) were unpaired; of these:
    9273865 (74.49%) aligned 0 times
    1058049 (8.50%) aligned exactly 1 time
    2117945 (17.01%) aligned >1 times
25.51% overall alignment rate

12006668 reads; of these:
  12006668 (100.00%) were unpaired; of these:
    8140818 (67.80%) aligned 0 times
    1016314 (8.46%) aligned exactly 1 time
    2849536 (23.73%) aligned >1 times
32.20% overall alignment rate

4539491 reads; of these:
  4539491 (100.00%) were unpaired; of these:
    1541573 (33.96%) aligned 0 times
    1398978 (30.82%) aligned exactly 1 time
    1598940 (35.22%) aligned >1 times
66.04% overall alignment rate

18303254 reads; of these:
  18303254 (100.00%) were unpaired; of these:
    11218117 (61.29%) aligned 0 times
    2862256 (15.64%) aligned exactly 1 time
    4222881 (23.07%) aligned >1 times
38.71% overall alignment rate

2443369 reads; of these:
  2443369 (100.00%) were unpaired; of these:
    960615 (39.32%) aligned 0 times
    396696 (16.24%) aligned exactly 1 time
    1086058 (44.45%) aligned >1 times
60.68% overall alignment rate

12065549 reads; of these:
  12065549 (100.00%) were unpaired; of these:
    9236644 (76.55%) aligned 0 times
    1113706 (9.23%) aligned exactly 1 time
    1715199 (14.22%) aligned >1 times
23.45% overall alignment rate
...
```


--------------------------------------------------------------------------------------------
## Need seq2iso tables for host and symbiont

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

#### de novo assemblies Siderastrea host, symbiont
Need hits seq name list (suffixed) for de novo assemblies of host and symbiont

cp Sid_Host-Sym_seq2seq_orig.txt Sid_Host-Sym_seq2seq_orig.tab

#### Breviolum assembly
File:  Copia de Symbio_Sider.fna.emapper.annotations

Extract suffixed Trinity names for now.
```
enable_lmod
module load container_env perl
grep -P -o "(\w+_\w+)" 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta | head
grep -P -o "(\w+_\w+)" 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta > 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed_HITS.txt
```
Sid_radians_BreviolumB5_Avila-Medina_Trinity-seq2seq.txt

cp Sid_radians_BreviolumB5_Avila-Medina_Trinity-seq2seq.txt Sid_radians_BreviolumB5_Avila-Medina_Trinity-seq2seq.tab

#### Cladocopium assembly
cp davies_cladeC_seq2iso_suffixed.txt davies_cladeC_seq2iso_suffixed.tab


-------------------------------------------------------------------------------------------
## Concatenate host and symbiont seq2 iso table

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/

*This needs to be updated to seq2iso table when we make it*
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Sid_Host-Sym_seq2seq_orig.tab

*This needs to be updated to seq2iso table*
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Sid_radians_BreviolumB5_Avila-Medina_Trinity-seq2seq.tab

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/davies_cladeC_seq2iso_suffixed.tab

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Sid_Host-Sym_seq2seq_orig.tab /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/Sid_radians_BreviolumB5_Avila-Medina_Trinity-seq2seq.tab /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/davies_cladeC_seq2iso_suffixed.tab > Sid_hybrid_seq2iso.tab
```
**For some reason the last line of file and first line of subsequent files were merged together?**
(I manually fixed this, but not sure why this happened)

--------------------------------------------------------------------------------------------
## Count expression 
Generate counts file for all samples

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/

```
nano countexpression_hybridref.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_hybridref.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_hybridref

enable_lmod
module load container_env python

python2 /cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/*_nof_hybridref.sam -g /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/Sid_hybrid_seq2iso.tab
```

```
sbatch countexpression_hybridref.sh
```


---------------------------------------------------------------------------------------------
## Parse expression to table

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref/

#### Make genelist .txt file
Format should be (one column):

GeneName |
--- |
gene1name |
gene2name |
gene3name |
gene4name |

```
nano ParseExpression.sh
```

```
#!/bin/bash -l

#SBATCH -o ParseExpression.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=ParseExpression

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/ParseExpression2BigTable_advbioinf.py Sid_genelist_seq-iso-combo.txt Sid_hybridref_counts.txt nomatch *_nof_hybridref_counts.txt
```

```
sbatch ParseExpression.sh
```

#### Output 
Sid_hybridref_counts.txt

--------------------------------------------------------------------------------------------
*Updates to symbiont references for hybrid ref mapping*

## Filter symbiont transcriptome references in same way as our de novo references

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/


### Clade B reference
The Avila-Medina reference has Trinity 'gene' names but is missing len= and path=.

head 173968_Breviolum_B5_Sid_radians_Avila-Medina.fna
```
>TRINITY_DN100007_c0_g1_i1
GTGATGTCAATGTCCTTTTGCACACGGATGATTTTGATAGGGTGTACCACTGGAAAATCCATCTGCATTT
GATATCCTAGTTCATCTTCAATTTCTTCATGTGGTTGAATCTGAAAGAATTGGAACAACCGATCCCCCTT
CGATCCGAAAATCCTGGTGCATTTGCCATCCATGCAAACATCAACGAACACCTCATTGCCAACAGACACA
>TRINITY_DN100007_c0_g2_i1
GTGATGTCAATGTCCTTTTGCACACGGATGATTTTGATAGGGTGTACCACTGGAAAATCCATCTGCATTT
GATATCCTAGTTCATCTTCAATTTCTTCATGTGATTGAATCTGAAAGAACTGGAACAACCGATCCCCCTT
CGATCCGAAAATCCTGGTGCATTTGCCATCCATGCAAACATCAACGAACACCTCATTGCCAACAGACACA
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta
```

```
The total number of sequences is 173968
The average sequence length is 627
The total number of bases is 109245170
The minimum sequence length is 201
The maximum sequence length is 2459
The N50 is 1040
Median Length = 258
contigs < 150bp = 0
contigs >= 500bp = 66248
contigs >= 1000bp = 36581
contigs >= 2000bp = 5807
```

#### Filter by a length threshold (> 500 bp)
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 _500lnThresh 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta
```

```
Number of total seqs for 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta: 173968
Number of seqs over 500 for 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta: 66249
```

**Best practice to rename with contig number**

```
grep -c ">" 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed_500lnThresh.fasta
```
66248

```
cp 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed_500lnThresh.fasta 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh.fasta
```

cp 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh.fasta 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh.txt


Import into R to obtain & add back in the length information into seq name


cp 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh_len-space.txt 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh_len.fasta


#### Longest Isoform

*Need length information in line with Trinity 'gene' name for the script below*
Standard Trinity.fasta output includes Trinity 'gene' name, length, and path (tab delimited)

**CONTINUE HERE WITH FASTA FILE WITH len=**

```
nano LongestIsoform.sh
```

```
#!/bin/bash -l

#SBATCH -o LongestIsoform.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=LongestIsoform

/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_longest_isoform_per_gene_trinity.py -i 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh_len.fasta -o Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
```

```
sbatch LongestIsoform.sh
```

grep -c ">" Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
66248

***^ DID THIS WORK OR NOT???***

Open txt file and use regular expressions:

Remove:  >TRINITY_

Find
(\w+_c\d_g\d)\w+\s\w+.\w+

Replace
$1

*13883 duplicate values found of Trinity 'genes'*
From 2-28 duplicates per given name_contig_gene
For example:  DN100580_c0_g1 showed up twice, so two isoforms existed for the same contig_gene

*Should get down to 52365 unique contig_genes (according to Excel remove duplicates)*

Maybe the _Symbiont suffix that I added is interfering with the script recognizing the isoforms? (Yes)

cp Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_noSuffix.txt 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_noSuffix.fasta

```
#!/bin/bash -l

#SBATCH -o LongestIsoform.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=LongestIsoform

/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_longest_isoform_per_gene_trinity.py -i 66248_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_noSuffix.fasta -o Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
```

*Double >> was introduced again*
```
>>TRINITY_DN714416_c0_g1_i1 len=1537
```

```
grep -c ">" Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
```
52508

grep -c "TRINITY" Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
52508

**Best practice to rename with contig number**

cp Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta

#### Check the >500 bp contig longest isoform assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta
```

```
The total number of sequences is 52508
The average sequence length is 1178
The total number of bases is 61857756
The minimum sequence length is 500
The maximum sequence length is 2459
The N50 is 1372
Median Length = 687
contigs < 150bp = 0
contigs >= 500bp = 52508
contigs >= 1000bp = 28949
contigs >= 2000bp = 4806
```

#### Edit name - add unique identifier for this specific reference

File:  52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta

#### Copy Trinity 'gene' contig names

*Need whole Trinity 'gene' contig name line with length, path (rather than contig name only)*
```
enable_lmod
module load container_env perl
grep -P -o "(\w+)\s\w+.\w+" 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta | head
grep -P -o "(\w+)\s\w+.\w+" 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta > 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_names.txt
```

In VS Studio regex
Find      ..TRINITY.(\w+)\s\w+.\w+
Replace   $1

### Rename contig names to something useful
Use script fasta_name_replacer.py

**Re-add suffix _Symbiont**

*Keep record of new names and corresponding original names*
- Create Tab delimited, 2 column table of OriginalName\tNewName (no headers).

head Breviolum_Avila-Medina_renamer.txt
```
>TRINITY_DN714416_c0_g1_i1 len=1537	Breviolum_Sid-radians_Avila-Medina_DN714416_c0_g1_i1_Symbiont
>TRINITY_DN167165_c0_g1_i1 len=516	Breviolum_Sid-radians_Avila-Medina_DN167165_c0_g1_i1_Symbiont
>TRINITY_DN767352_c0_g1_i1 len=1890	Breviolum_Sid-radians_Avila-Medina_DN767352_c0_g1_i1_Symbiont
>TRINITY_DN58758_c0_g2_i1 len=1296	Breviolum_Sid-radians_Avila-Medina_DN58758_c0_g2_i1_Symbiont
>TRINITY_DN128483_c0_g1_i1 len=543	Breviolum_Sid-radians_Avila-Medina_DN128483_c0_g1_i1_Symbiont
>TRINITY_DN614257_c0_g1_i1 len=644	Breviolum_Sid-radians_Avila-Medina_DN614257_c0_g1_i1_Symbiont
>TRINITY_DN258489_c0_g1_i1 len=551	Breviolum_Sid-radians_Avila-Medina_DN258489_c0_g1_i1_Symbiont
>TRINITY_DN7697_c0_g2_i1 len=1243	Breviolum_Sid-radians_Avila-Medina_DN7697_c0_g2_i1_Symbiont
>TRINITY_DN435435_c1_g4_i2 len=730	Breviolum_Sid-radians_Avila-Medina_DN435435_c1_g4_i2_Symbiont
>TRINITY_DN847198_c0_g2_i1 len=556	Breviolum_Sid-radians_Avila-Medina_DN847198_c0_g2_i1_Symbiont
```

```
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/fasta_name_replacer.py -i 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso.fasta -n Breviolum_Avila-Medina_renamer.txt -o 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta
```

head 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta
```
>Breviolum_Sid-radians_Avila-Medina_DN714416_c0_g1_i1_Symbiont
GGCACAAGTAAAAAAAATCTTCGAAAGCAAATACAAGACGATGAAAGAAGAGCAGTTTAGTATATTTTGGTAAACTACTGTAACTTAAACTAAACAGGGTGTCGATTAGCACTTGCACTAAATACCACCGGTTTTTTTTAGGTGGTATCAGACAATCATAAAAATGTAAAATTTTTTTGTACTA
```

*The file 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta can now be used to make the hybrid reference and mapping.*



### Clade C reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

Specific for Davies reference because different contig name format (compared to standard Trinity output).
An approach for subsetting the longest contig per gene from the Davies transcriptome using her seq2iso.tab file and the input .fasta file.

```
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/fasta_longest_contig_per_gene_seq2genetable.py -i 65838_davies_cladeC_feb_fixed-clean_suffixed.fasta -s davies_cladeC_seq2iso_suffixed.tab -o davies_cladeC_LongestContig.fasta
```

#### Filter by a length threshold (> 500 bp)

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 500lnThresh davies_cladeC_LongestContig.fasta
```

```
Number of total seqs for davies_cladeC_LongestContig.fasta: 47264
Number of seqs over 500 for davies_cladeC_LongestContig.fasta: 47263
```

```
grep -c '>' davies_cladeC_LongestContig500lnThresh.fasta
```
47262

**Best practice to rename with contig number**
```
cp davies_cladeC_LongestContig500lnThresh.fasta 47262_davies_cladeC_LongestContig_500lnThresh.fasta
```

#### Check assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 47262_davies_cladeC_LongestContig_500lnThresh.fasta
```

```
The total number of sequences is 47262
The average sequence length is 1364
The total number of bases is 64485324
The minimum sequence length is 500
The maximum sequence length is 18168
The N50 is 1588
Median Length = 1519
contigs < 150bp = 0
contigs >= 500bp = 47262
contigs >= 1000bp = 27124
contigs >= 2000bp = 7623
```

*Also can be done in R*
```
GeneTable <- read.delim(file="/Users/danbarshis/dansstuff/Projeks/ODU/Projeks/reference_omes/CladeC_Symbiodinium_transcriptome/daviesSeq2IsoPlusLengthandsorted.txt")

GeneTable <- GeneTable[order(GeneTable$GeneName, GeneTable$Length, decreasing=c(FALSE,TRUE)),]

head(!(duplicated(GeneTable$GeneName)))

write.table(GeneTable[!(duplicated(GeneTable$GeneName)),],file="/Users/danbarshis/dansstuff/Projeks/ODU/Projeks/reference_omes/CladeC_Symbiodinium_transcriptome/LongestContigperGene.txt", sep='\t', row.names=F, quote=F)
```

#### Edit name - add unique identifier for each specific reference

File:  47262_davies_cladeC_LongestContig_500lnThresh.fasta

>comp105411_c0_seq1_Symbiont

#### Copy Trinity 'gene' contig names

*Need whole Trinity 'gene' contig name line with length, path (rather than contig name only)*
```
enable_lmod
module load container_env perl
grep -P -o "(\w+.Symbiont)" 47262_davies_cladeC_LongestContig_500lnThresh.fasta | head
grep -P -o "(\w+.Symbiont)" 47262_davies_cladeC_LongestContig_500lnThresh.fasta > 47262_davies_cladeC_LongestContig_500lnThresh_names.txt
```

### Rename contig names to something useful
Use script fasta_name_replacer.py

*Keep record of new names and corresponding original names*
- Create Tab delimited, 2 column table of OriginalName\tNewName (no headers).

head Cladocopium_Davies_renamer.txt
```
comp105411_c0_seq1_Symbiont	Cladocopium_Davies_comp105411_c0_seq1_Symbiont
comp72706_c0_seq1_Symbiont	Cladocopium_Davies_comp72706_c0_seq1_Symbiont
comp32947_c0_seq1_Symbiont	Cladocopium_Davies_comp32947_c0_seq1_Symbiont
comp187852_c0_seq1_Symbiont	Cladocopium_Davies_comp187852_c0_seq1_Symbiont
comp298882_c0_seq1_Symbiont	Cladocopium_Davies_comp298882_c0_seq1_Symbiont
comp460517_c0_seq1_Symbiont	Cladocopium_Davies_comp460517_c0_seq1_Symbiont
comp110971_c0_seq2_Symbiont	Cladocopium_Davies_comp110971_c0_seq2_Symbiont
comp268310_c0_seq1_Symbiont	Cladocopium_Davies_comp268310_c0_seq1_Symbiont
comp205538_c0_seq1_Symbiont	Cladocopium_Davies_comp205538_c0_seq1_Symbiont
comp169718_c0_seq1_Symbiont	Cladocopium_Davies_comp169718_c0_seq1_Symbiont
```

```
python2 /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/fasta_name_replacer.py -i 47262_davies_cladeC_LongestContig_500lnThresh.fasta -n Cladocopium_Davies_renamer.txt -o 47262_Davies_Cladocopium_LongestContig-500_renamed.fasta
```

head 47262_Davies_Cladocopium_LongestContig-500_renamed.fasta
```
>Cladocopium_Davies_comp105411_c0_seq1_Symbiont
GGATGTGTGCAACAACTTCCTGGAGTCGCATCATATTACCCCGGAGCTCTCCATCCGCATGAAGAAGTTCGTAGTGTCGTGTCACCAGCAATCCTGGCTGGAGAGCAAGCTCAAGGAGGAGCAGCAGCTTTTGTCGAAGATGCCGGAGATCTTGCAGTCAGATTTGTACGAGGAGTCCAGGGGCCCGAACCTCTGTAAGAGCACCTTCTGGAATGAGTTTCAGGAGAAGTTCCGCGCCTGTTTCCGCCACGTTTGCATGGACGCGGTGACGGAGGTCGTGGCACAAAAGCACGACATGATCTTCACCGATTTCCATCACGCGCATTCCACGCTCTACATTTTATCGGGCTCGTTCTGCTACACGCTGGCGAGGAAGCACACCTTGTTGGGGCAGCTGGGCCAACTGCTGGCTCCCCGACATGTGCGCTCACGCAAGGTGAAGGTGGAAGGTGGTCGGAGCATCTGTGAGATTGCCCTCTGGACCAACTGGATCCACACAGGCTCCCTACAGGTGCTGAAGAGCGGCACCTACCTGTCGGTGAACGTGGAAGTGGCCAGCGACATCATCGCTGCCTATCCGGAAGCCCAGCGCACCGCGCACGCGTACGGTCGCATGGTGGTGATCCAGCTCCATCACTGCAAGAACGATTGCACGGACCTGACGCCTTTGGAGATAGACTTCAAGTCGCTGAGGAAGGTGACGGTGAAGTTCATGACAGAGGGCCACTGCATCTTCCTGACGCACTTCAAAAAGGAGGCAGGGACAGAGGCGGCGCTCATGCAGGAAGCCCTTACGAGCAAGATTCAGTCAGATCCACAAAATCCGGCACACTACTTGGAGCACCCGGTCTTCTTGGACACTGAGAATCTGGAAGACCTAGCCAAGTTGAGGGACCACATCCTGGCCAGCCGGAACTTGGTGGTCCTGTTGACACCCGGCATCTTCGAACGGCCGTGGTGCTTGGTGGAGATCGTCACGGCCTTCAAGCAAAACCGGAACATTGTGTTAGTTGAAATCCAAAAGCCGGACATGAAGTTCGAGTATCCGGATGAGGAGTTTATCCTTGACCTTTGTGACGGGAGACTCTTCAGCGAAACGGATATGCAAGTCTTTGCACTTGAGCAGATTGAGTTGACGGACATTGAGGTGGCCTTGCGACACGTGCTCACAAAGATCGCCCTGCCGTTTTCGCCACACAAAGCGAG
>Cladocopium_Davies_comp72706_c0_seq1_Symbiont
TTTTCACCGAGCTCTTCAAGCTGGCCTTAAGCTTCGCGTTCCTATCGCTGGAGA
```

*The file 47262_Davies_Cladocopium_LongestContig-500_renamed.fasta can now be used to make the hybrid reference and mapping.*


--------------------------------------------------------------------------------------------
## Create hybrid reference for mapping

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/

### Symbiont
52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta

47262_Davies_Cladocopium_LongestContig-500_renamed.fasta

### Host
*De novo assembly*
19222_Sid_GoodCoral_500lnThresh_Final.fasta

#### Copy Host de novo assembly & symbiont references
```
cp /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/19222_Sid_GoodCoral_500lnThresh_Final.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/47262_Davies_Cladocopium_LongestContig-500_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/
```

--------------------------------------------------------------------------------------------
## Concatenate host reference and symbiont references to make mapping file

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/19222_Sid_GoodCoral_500lnThresh_Final.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/47262_Davies_Cladocopium_LongestContig-500_renamed.fasta > Sid_hybridref_final.fasta
```

## Check assembly details

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Sid_hybridref_final.fasta
```

```
The total number of sequences is 118992
The average sequence length is 1243
The total number of bases is 148010700
The minimum sequence length is 500
The maximum sequence length is 32630
The N50 is 1439
Median Length = 1570
contigs < 150bp = 0
contigs >= 500bp = 118992
contigs >= 1000bp = 62766
contigs >= 2000bp = 14258
```

--------------------------------------------------------------------------------------------
## Map all the samples to concatenated (host plus symbiont) hybrid reference transcriptome

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/Sid_hybridref_final.fasta hybridref_final
```

#### Mapping to hybridreference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref_final/

```
nano mapreads_hybridref_final.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_hybridref_final.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_hybridref

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x hybridref_final -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_hybridref.sam -k 5\n; done
```

```
sbatch mapreads_hybridref_final.sh
```

```
cat bowtie2_hybridref_final.txt
```

```
10127834 reads; of these:
  10127834 (100.00%) were unpaired; of these:
    6789363 (67.04%) aligned 0 times
    1937909 (19.13%) aligned exactly 1 time
    1400562 (13.83%) aligned >1 times
32.96% overall alignment rate
4535297 reads; of these:
  4535297 (100.00%) were unpaired; of these:
    2702008 (59.58%) aligned 0 times
    904330 (19.94%) aligned exactly 1 time
    928959 (20.48%) aligned >1 times
40.42% overall alignment rate
12959491 reads; of these:
  12959491 (100.00%) were unpaired; of these:
    9222597 (71.16%) aligned 0 times
    2408083 (18.58%) aligned exactly 1 time
    1328811 (10.25%) aligned >1 times
28.84% overall alignment rate
14417329 reads; of these:
  14417329 (100.00%) were unpaired; of these:
    9857318 (68.37%) aligned 0 times
    2837642 (19.68%) aligned exactly 1 time
    1722369 (11.95%) aligned >1 times
31.63% overall alignment rate
15619667 reads; of these:
  15619667 (100.00%) were unpaired; of these:
    10652696 (68.20%) aligned 0 times
    3041671 (19.47%) aligned exactly 1 time
    1925300 (12.33%) aligned >1 times
31.80% overall alignment rate
22723159 reads; of these:
  22723159 (100.00%) were unpaired; of these:
    15190421 (66.85%) aligned 0 times
    4622015 (20.34%) aligned exactly 1 time
    2910723 (12.81%) aligned >1 times
33.15% overall alignment rate
23167285 reads; of these:
  23167285 (100.00%) were unpaired; of these:
    16150939 (69.71%) aligned 0 times
    4497853 (19.41%) aligned exactly 1 time
    2518493 (10.87%) aligned >1 times
30.29% overall alignment rate
19306925 reads; of these:
  19306925 (100.00%) were unpaired; of these:
    12771136 (66.15%) aligned 0 times
    4190751 (21.71%) aligned exactly 1 time
    2345038 (12.15%) aligned >1 times
33.85% overall alignment rate
21243157 reads; of these:
  21243157 (100.00%) were unpaired; of these:
    14487991 (68.20%) aligned 0 times
    4115647 (19.37%) aligned exactly 1 time
    2639519 (12.43%) aligned >1 times
31.80% overall alignment rate
20103321 reads; of these:
  20103321 (100.00%) were unpaired; of these:
    16551673 (82.33%) aligned 0 times
    2441154 (12.14%) aligned exactly 1 time
    1110494 (5.52%) aligned >1 times
17.67% overall alignment rate
18367528 reads; of these:
  18367528 (100.00%) were unpaired; of these:
    12091605 (65.83%) aligned 0 times
    3832069 (20.86%) aligned exactly 1 time
    2443854 (13.31%) aligned >1 times
34.17% overall alignment rate
22832214 reads; of these:
  22832214 (100.00%) were unpaired; of these:
    14968323 (65.56%) aligned 0 times
    4851909 (21.25%) aligned exactly 1 time
    3011982 (13.19%) aligned >1 times
34.44% overall alignment rate
12449859 reads; of these:
  12449859 (100.00%) were unpaired; of these:
    10278912 (82.56%) aligned 0 times
    1185367 (9.52%) aligned exactly 1 time
    985580 (7.92%) aligned >1 times
17.44% overall alignment rate
12316099 reads; of these:
  12316099 (100.00%) were unpaired; of these:
    8247259 (66.96%) aligned 0 times
    2374825 (19.28%) aligned exactly 1 time
    1694015 (13.75%) aligned >1 times
33.04% overall alignment rate
21586925 reads; of these:
  21586925 (100.00%) were unpaired; of these:
    14577368 (67.53%) aligned 0 times
    4294950 (19.90%) aligned exactly 1 time
    2714607 (12.58%) aligned >1 times
32.47% overall alignment rate
12103752 reads; of these:
  12103752 (100.00%) were unpaired; of these:
    9737342 (80.45%) aligned 0 times
    1318030 (10.89%) aligned exactly 1 time
    1048380 (8.66%) aligned >1 times
19.55% overall alignment rate
8479518 reads; of these:
  8479518 (100.00%) were unpaired; of these:
    6338762 (74.75%) aligned 0 times
    1384531 (16.33%) aligned exactly 1 time
    756225 (8.92%) aligned >1 times
25.25% overall alignment rate
21018098 reads; of these:
  21018098 (100.00%) were unpaired; of these:
    14237120 (67.74%) aligned 0 times
    4136326 (19.68%) aligned exactly 1 time
    2644652 (12.58%) aligned >1 times
32.26% overall alignment rate
7890229 reads; of these:
  7890229 (100.00%) were unpaired; of these:
    5095589 (64.58%) aligned 0 times
    1300507 (16.48%) aligned exactly 1 time
    1494133 (18.94%) aligned >1 times
35.42% overall alignment rate
18032416 reads; of these:
  18032416 (100.00%) were unpaired; of these:
    12595702 (69.85%) aligned 0 times
    3486649 (19.34%) aligned exactly 1 time
    1950065 (10.81%) aligned >1 times
30.15% overall alignment rate
8523378 reads; of these:
  8523378 (100.00%) were unpaired; of these:
    5633695 (66.10%) aligned 0 times
    1617779 (18.98%) aligned exactly 1 time
    1271904 (14.92%) aligned >1 times
33.90% overall alignment rate
12006668 reads; of these:
  12006668 (100.00%) were unpaired; of these:
    9045482 (75.34%) aligned 0 times
    1151596 (9.59%) aligned exactly 1 time
    1809590 (15.07%) aligned >1 times
24.66% overall alignment rate
19967419 reads; of these:
  19967419 (100.00%) were unpaired; of these:
    13690050 (68.56%) aligned 0 times
    3872228 (19.39%) aligned exactly 1 time
    2405141 (12.05%) aligned >1 times
31.44% overall alignment rate
4539491 reads; of these:
  4539491 (100.00%) were unpaired; of these:
    2075048 (45.71%) aligned 0 times
    1379143 (30.38%) aligned exactly 1 time
    1085300 (23.91%) aligned >1 times
54.29% overall alignment rate
10328100 reads; of these:
  10328100 (100.00%) were unpaired; of these:
    6427057 (62.23%) aligned 0 times
    1930840 (18.70%) aligned exactly 1 time
    1970203 (19.08%) aligned >1 times
37.77% overall alignment rate
11972400 reads; of these:
  11972400 (100.00%) were unpaired; of these:
    8175152 (68.28%) aligned 0 times
    2311497 (19.31%) aligned exactly 1 time
    1485751 (12.41%) aligned >1 times
31.72% overall alignment rate
17479996 reads; of these:
  17479996 (100.00%) were unpaired; of these:
    12000461 (68.65%) aligned 0 times
    3519965 (20.14%) aligned exactly 1 time
    1959570 (11.21%) aligned >1 times
31.35% overall alignment rate
15309393 reads; of these:
  15309393 (100.00%) were unpaired; of these:
    9859042 (64.40%) aligned 0 times
    2906780 (18.99%) aligned exactly 1 time
    2543571 (16.61%) aligned >1 times
35.60% overall alignment rate
17272901 reads; of these:
  17272901 (100.00%) were unpaired; of these:
    11628569 (67.32%) aligned 0 times
    3440037 (19.92%) aligned exactly 1 time
    2204295 (12.76%) aligned >1 times
32.68% overall alignment rate
14020924 reads; of these:
  14020924 (100.00%) were unpaired; of these:
    9592761 (68.42%) aligned 0 times
    2625660 (18.73%) aligned exactly 1 time
    1802503 (12.86%) aligned >1 times
31.58% overall alignment rate
16213720 reads; of these:
  16213720 (100.00%) were unpaired; of these:
    11168389 (68.88%) aligned 0 times
    2979747 (18.38%) aligned exactly 1 time
    2065584 (12.74%) aligned >1 times
31.12% overall alignment rate
25535025 reads; of these:
  25535025 (100.00%) were unpaired; of these:
    17184162 (67.30%) aligned 0 times
    5029848 (19.70%) aligned exactly 1 time
    3321015 (13.01%) aligned >1 times
32.70% overall alignment rate
23328804 reads; of these:
  23328804 (100.00%) were unpaired; of these:
    15783362 (67.66%) aligned 0 times
    4564818 (19.57%) aligned exactly 1 time
    2980624 (12.78%) aligned >1 times
32.34% overall alignment rate
13100851 reads; of these:
  13100851 (100.00%) were unpaired; of these:
    8887387 (67.84%) aligned 0 times
    2547892 (19.45%) aligned exactly 1 time
    1665572 (12.71%) aligned >1 times
32.16% overall alignment rate
20666252 reads; of these:
  20666252 (100.00%) were unpaired; of these:
    14095990 (68.21%) aligned 0 times
    4243106 (20.53%) aligned exactly 1 time
    2327156 (11.26%) aligned >1 times
31.79% overall alignment rate
18750534 reads; of these:
  18750534 (100.00%) were unpaired; of these:
    12686938 (67.66%) aligned 0 times
    3865553 (20.62%) aligned exactly 1 time
    2198043 (11.72%) aligned >1 times
32.34% overall alignment rate
19423159 reads; of these:
  19423159 (100.00%) were unpaired; of these:
    13254138 (68.24%) aligned 0 times
    4061712 (20.91%) aligned exactly 1 time
    2107309 (10.85%) aligned >1 times
31.76% overall alignment rate
19253832 reads; of these:
  19253832 (100.00%) were unpaired; of these:
    13258837 (68.86%) aligned 0 times
    3759940 (19.53%) aligned exactly 1 time
    2235055 (11.61%) aligned >1 times
31.14% overall alignment rate
16886305 reads; of these:
  16886305 (100.00%) were unpaired; of these:
    11450734 (67.81%) aligned 0 times
    3356240 (19.88%) aligned exactly 1 time
    2079331 (12.31%) aligned >1 times
32.19% overall alignment rate
17696748 reads; of these:
  17696748 (100.00%) were unpaired; of these:
    12330955 (69.68%) aligned 0 times
    3638606 (20.56%) aligned exactly 1 time
    1727187 (9.76%) aligned >1 times
30.32% overall alignment rate
20281418 reads; of these:
  20281418 (100.00%) were unpaired; of these:
    14470130 (71.35%) aligned 0 times
    3858409 (19.02%) aligned exactly 1 time
    1952879 (9.63%) aligned >1 times
28.65% overall alignment rate
18064774 reads; of these:
  18064774 (100.00%) were unpaired; of these:
    15050218 (83.31%) aligned 0 times
    1817736 (10.06%) aligned exactly 1 time
    1196820 (6.63%) aligned >1 times
16.69% overall alignment rate
18881427 reads; of these:
  18881427 (100.00%) were unpaired; of these:
    12755172 (67.55%) aligned 0 times
    3783183 (20.04%) aligned exactly 1 time
    2343072 (12.41%) aligned >1 times
32.45% overall alignment rate
18303254 reads; of these:
  18303254 (100.00%) were unpaired; of these:
    13764456 (75.20%) aligned 0 times
    2766631 (15.12%) aligned exactly 1 time
    1772167 (9.68%) aligned >1 times
24.80% overall alignment rate
2443369 reads; of these:
  2443369 (100.00%) were unpaired; of these:
    1193725 (48.86%) aligned 0 times
    417655 (17.09%) aligned exactly 1 time
    831989 (34.05%) aligned >1 times
51.14% overall alignment rate
23501378 reads; of these:
  23501378 (100.00%) were unpaired; of these:
    16062503 (68.35%) aligned 0 times
    4817294 (20.50%) aligned exactly 1 time
    2621581 (11.16%) aligned >1 times
31.65% overall alignment rate
13704567 reads; of these:
  13704567 (100.00%) were unpaired; of these:
    8963366 (65.40%) aligned 0 times
    3042827 (22.20%) aligned exactly 1 time
    1698374 (12.39%) aligned >1 times
34.60% overall alignment rate
13293412 reads; of these:
  13293412 (100.00%) were unpaired; of these:
    9066711 (68.20%) aligned 0 times
    2779065 (20.91%) aligned exactly 1 time
    1447636 (10.89%) aligned >1 times
31.80% overall alignment rate
14496181 reads; of these:
  14496181 (100.00%) were unpaired; of these:
    10040683 (69.26%) aligned 0 times
    2965436 (20.46%) aligned exactly 1 time
    1490062 (10.28%) aligned >1 times
30.74% overall alignment rate
16287244 reads; of these:
  16287244 (100.00%) were unpaired; of these:
    10913795 (67.01%) aligned 0 times
    3587925 (22.03%) aligned exactly 1 time
    1785524 (10.96%) aligned >1 times
32.99% overall alignment rate
16428742 reads; of these:
  16428742 (100.00%) were unpaired; of these:
    12139751 (73.89%) aligned 0 times
    2743989 (16.70%) aligned exactly 1 time
    1545002 (9.40%) aligned >1 times
26.11% overall alignment rate
12065549 reads; of these:
  12065549 (100.00%) were unpaired; of these:
    10159492 (84.20%) aligned 0 times
    1260757 (10.45%) aligned exactly 1 time
    645300 (5.35%) aligned >1 times
15.80% overall alignment rate
18031756 reads; of these:
  18031756 (100.00%) were unpaired; of these:
    12239898 (67.88%) aligned 0 times
    3399278 (18.85%) aligned exactly 1 time
    2392580 (13.27%) aligned >1 times
32.12% overall alignment rate
29095469 reads; of these:
  29095469 (100.00%) were unpaired; of these:
    24065590 (82.71%) aligned 0 times
    3140200 (10.79%) aligned exactly 1 time
    1889679 (6.49%) aligned >1 times
17.29% overall alignment rate
21883993 reads; of these:
  21883993 (100.00%) were unpaired; of these:
    15176420 (69.35%) aligned 0 times
    4248851 (19.42%) aligned exactly 1 time
    2458722 (11.24%) aligned >1 times
30.65% overall alignment rate
23679838 reads; of these:
  23679838 (100.00%) were unpaired; of these:
    17311958 (73.11%) aligned 0 times
    3997153 (16.88%) aligned exactly 1 time
    2370727 (10.01%) aligned >1 times
26.89% overall alignment rate
```


--------------------------------------------------------------------------------------------
## Count expression 
Generate counts file for all samples

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref_final/

*Do not need seq2iso table because it is contig based*

```
nano countexpression_hybridref.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_hybridref.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_hybridref

enable_lmod
module load container_env python

python2 /cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref_final/*_nof_hybridref.sam
```

```
sbatch countexpression_hybridref.sh
```


---------------------------------------------------------------------------------------------
## Parse expression to table

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/mapped_hybridref_final/

*In this case it is contig name list (in the order that the hybridref was concatenated*
19222_Sid_GoodCoral_500lnThresh_Final.fasta
52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta
47262_Davies_Cladocopium_LongestContig-500_renamed.fasta

```
grep -o -E "(Siderastrea_\w+.\w+)" 19222_Sid_GoodCoral_500lnThresh_Final.fasta | head
grep -o -E "(Siderastrea_\w+.\w+)" 19222_Sid_GoodCoral_500lnThresh_Final.fasta > 19222_Sid_GoodCoral_500lnThresh_Final_contigNames.txt
```

File:  ContigList_Sid_hybridref.txt

#### Make contig list .txt file
Format should be (one column):

GeneName |
--- |
gene1name |
gene2name |
gene3name |
gene4name |

```
nano ParseExpression.sh
```

```
#!/bin/bash -l

#SBATCH -o ParseExpression.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=ParseExpression

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/ParseExpression2BigTable_advbioinf.py ContigList_Sid_hybridref.txt Sid_hybridref_final_counts.txt nomatch *_nof_hybridref_counts.txt
```

```
sbatch ParseExpression.sh
```


#### Output 
Sid_hybridref_final_counts.txt

--------------------------------------------------------------------------------------------

no match
Breviolum_Sid-radians_Avila-Medina_DN100009_c1_g1_i1_Symbiont

head -4 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta
Breviolum_Sid-radians_Avila-Medina_DN714416_c0_g1_i1_Symbiont

