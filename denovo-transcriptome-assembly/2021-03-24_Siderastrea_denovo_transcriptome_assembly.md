---
layout: post
title: Siderastrea de novo transcriptome assembly
date: '2021-01-06'
categories: Protocols
tags: [RNASeq, Bioinformatics]
---


# *Siderastrea* de novo transcriptome assembly
## Build assembly

### *Siderastrea siderea*
- Puerto Morelos, Mexico
- Samples were collected at two time points, Year 1 and Year 2 (of transplant)
- Ana Martinez (PhD candidate) & Adina Paytan (PI) & Dan Barshis (PI)

### Analyzed by Veronica Radice, Barshis Lab, Old Dominion University

## Data (Year 1 samples)
- Paired end 150bp lane (12 Sid samples)  
- Single end 150bp lane (17 Sid samples)

### Files:
Raw data (Year 1 samples)
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/raw_Year1/

Clipped, trimmed, nofilter files (Paired, R1 and R2; and single end, R1 only):
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/Year1_clip-trim_nofilter/

Obtained from:
> /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/

*Use the stillpaired files here as the paired files (folder also has singletons):*
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/Year1_clip-trim_nofilter_paired/

Obtained from: 
> /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/


### Raw data backed up on RC drive:
Contains data for multiple species:  *Porites astreoides*, *Porites porites*, *Siderastrea siderea*)

Siderastrea, **Single end 50bp** (Year 2) (for gene expression analysis):
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2017_April/gslserver.qb3.berkeley.edu/170419_50SR_HS4K2A/Paytan/
Processed:
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/Siderastrea/clippedtrimmed_nofilter/

Siderastrea, **Single end 150bp** (n=11) (Year 1): 
APDB13
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/141013_HS3B/Project_Paytan/

Within one given folder:
- APDB13A_index2_CGATGT_L002_R1_001.fastq.gz
- APDB13A_index2_CGATGT_L002_R1_002.fastq.gz
- APDB13A_index2_CGATGT_L002_R1_003.fastq.gz
- APDB13A_index2_CGATGT_L002_R1_004.fastq.gz


Siderastrea, **Paired end 150bp** (n=12) (Year 1):
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/141015_HS2B/Project_Paytan/

Porites porites, **Paired end 150bp** (n=11) (Year 1):
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/141015_HS2B/Project_Paytan/

Porites astreoides, **Single end 150bp** (n=11) (Year 1):
APBD14
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/141018_HS3B/Project_PAYTAN/

Mixed corals, **Single end 150bp** (n=11) (Year 1):
APBD15
> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/141018_HS3B/Project_PAYTAN/

#### Analyses originally run here
> /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/

#### Moved here
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/


----------------------------------------------------------------------------------------------
## Quality Control – Run QC on the FastQ files

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/fastqc/

fastqc/1.9 (Default)
```
enable_lmod
module load container_env trinity
module avail fastqc
module load fastqc/1.9
crun fastqc /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/*_Sid_yr2_R1_clippedtrimmed_nofilter.fastq
```

You can open this in a web browser. 
*_fastqc.html
The & on the end of the commend means that although firefox will launch, it won’t stop you from running other programs in your shell.
```
firefox PC_349_C_Sid_yr2_R1_clippedtrimmed_nofilter_fastqc.html &
```


----------------------------------------------------------------------------------------------
## Best Practices for De Novo Transcriptome Assembly with Trinity
- Despite a steady increase in the availablity of tools and documented pipelines for building transcriptome assemblies, de novo transcriptome assembly from relative short Illumina paired-end reads remains an extremely challenging endeavor.
- While early genome assemblers used pairwise overlaps between long reads to extend contigs, this approach is unfeasible when dealing with hundreds of millions of reads. 
- Thus, de novo transcriptome asssemblers use DeBruijn graphs, which are constructed and extended based upon kmers, i.e. subsequence of length k found in reads. 
- While this makes the assembly process computationally tractable, it can lead to fragmented assemblies of a large number of contigs that are subsequences of the underlying true transcripts. 
- Some of the factors that lead to this fragmentation are sequencing errors, polymorphism, sequence repeats, and for more lowly expressed transcripts, stochasticicty of read depth that leads to gaps in coverage. 
- Transcriptome assemblers, unlike genome assemblers, must handle the wide range of depth of coverage due to gene expression variation. 
- Our goal in developing a best practice pipeline is to produce most contiguous, error-free and complete transcriptome assemblies given these challenges. While a specific pipeline may produce near optimal results in most scenarios, there may be particular biological and technical factors worth considering that may lead to modest changes.

[https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html)


----------------------------------------------------------------------------------------------
### Use Trinity to create reference assembly

Short video overview of Trinity algorithm:
[https://www.broadinstitute.org/videos/trinity-how-it-works](https://www.broadinstitute.org/videos/trinity-how-it-works)

- Trinity assembles transcript sequences from Illumina RNA-Seq data
[https://github.com/trinityrnaseq/trinityrnaseq/wiki](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- Trinity, developed at the Broad Institute and the Hebrew University of Jerusalem, represents a novel method for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data 
- Trinity combines three independent software modules: Inchworm, Chrysalis, and Butterfly, applied sequentially to process large volumes of RNA-seq reads
- Trinity can be referenced as:
    Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. [doi.org/10.1038/nbt.1883](doi.org/10.1038/nbt.1883)
- Protocol for using Trinity for de novo transcriptome assembly and downstream analyses:
    Haas BJ, Papanicolaou A, Yassour M, Grabherr M, Blood PD, Bowden J, Couger MB, Eccles D, Li B, Lieber M, Macmanes MD, Ott M, Orvis J, Pochet N, Strozzi F, Weeks N, Westerman R, William T, Dewey CN, Henschel R, Leduc RD, Friedman N, Regev A. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nat Protoc. 2013 Aug;8(8):1494-512.
    [doi.org/10.1038/nprot.2013.084](doi.org/10.1038/nprot.2013.084)
- Performance tuning of Trinity is described in:
    Henschel R, Lieber M, Wu L, Nista, PM, Haas BJ, LeDuc R. Trinity RNA-Seq assembler performance optimization. XSEDE 2012 Proceedings of the 1st Conference of the Extreme Science and Engineering Discovery Environment: Bridging from the eXtreme to the campus and beyond. [doi.org/10.1145/2335755.2335842](doi.org/10.1145/2335755.2335842)


----------------------------------------------------------------------------------------------
Path to Trinity program:
> /opt/trinity/bin/

- Trinity includes additional options to automate various aspects of RNA-Seq read processing that should be considered prior to executing the de novo assembly
    - This includes quality trimming of reads using Trimmomatic, to reduce the number of reads that are subject to de novo assembly, improving on assembly run-time
    - NB:  I already have clipped and trimmed (no filter) fastq files ready to go
    - clippedtrimmed_nofilter.fastq files are output from the Trimclipfilterstatsbatch_advbioinf.py script


----------------------------------------------------------------------------------------------
#### Explaining the identifiers: Genes vs. Transcripts

File:  Trinity.fasta

- Trinity groups transcripts into clusters based on shared sequence content. 
- Such a transcript cluster is very loosely referred to as a 'gene'. 
- This information is encoded in the Trinity fasta accession.

Example:
```
>TRINITY_DN1000_c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
 AATCTTTTTTGGTATTGGCAGTACTGTGCTCTGGGTAGTGATTAGGGCAAAAGAAGACAC
 ACAATAAAGAACCAGGTGTTAGACGTCAGCAAGTCAAGGCCTTGGTTCTCAGCAGACAGA
 AGACAGCCCTTCTCAATCCTCATCCCTTCCCTGAACAGACATGTCTTCTGCAAGCTTCTC
 CAAGTCAGTTGTTCACAGGAACATCATCAGAATAAATTTGAAATTATGATTAGTATCTGA
 TAAAGCA
```

- The accession encodes the Trinity 'gene' and 'isoform' information. 
- In the example above, the accession 'TRINITY_DN1000_c115_g5_i1' indicates Trinity read cluster 'TRINITY_DN1000_c115', gene 'g5', and isoform 'i1'. 
- Because a given run of trinity involves many many clusters of reads, each of which are assembled separately, and because the 'gene' numberings are unique within a given processed read cluster, the 'gene' identifier should be considered an aggregate of the read cluster and corresponding gene identifier, which in this case would be 'TRINITY_DN1000_c115_g5'.
- So, in summary, the above example corresponds to 'gene id: TRINITY_DN1000_c115_g5' encoding 'isoform id: TRINITY_DN1000_c115_g5_i1'.


----------------------------------------------------------------------------------------------
## Control population assembly

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

The origin of all samples was a control site (Norte control or Pargos control):

NC_377A_C_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq
NC_378A_La_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq
NC_379A_N_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq
PC_373A_N_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq
PC_365A_C_Sid_yr1_R1_clippedtrimmed_nofilter.fastq
PC_366A_La_Sid_yr1_R1_clippedtrimmed_nofilter.fastq
PC_350A_La_Sid_yr1_R1_clippedtrimmed_nofilter.fastq
PC_364A_N_Sid_yr1_R1_clippedtrimmed_nofilter.fastq
PC_368A_C_Sid_yr1_R1_clippedtrimmed_nofilter.fastq

Assembly without orphaned singles (*_clippedtrimmed_nofilter_singles.fastq).

There are many R2 reads in the single (orphaned) files left over from the split of the files into R1 and R2 ("still_paired")
However, Trinity only accepts R1 sequences for single unpaired files (would need to add /1 suffix).

Run mixed PE and SE assembly without any "singles" (orphaned) files [leftover from "still_paired" files] because within these files there are some R2 sequences remaining.

For example:
```
grep "2:N:0:CCGTCC" NC_377A_C_Sid_yr1_R1_clippedtrimmed_nofilter_singles.fastq
```

Single (orphaned) files have much less reads than the paired read files:
```
grep -c @HS3:541 NC_377A_C_Sid_yr1_R1_clippedtrimmed_nofilter_singles.fastq
```
1,089,403

R1
```
grep -c @HS3:541 NC_377A_C_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq
```
8,738,919

R2
```
grep -c @HS3:541 NC_377A_C_Sid_yr1_R2_clippedtrimmed_nofilter_stillpaired.fastq
```
8,738,919

From Wiki:
If you have both paired and unpaired data, and the data are NOT strand-specific, you can combine the unpaired data with the left reads of the paired fragments. Be sure that the unpaired reads have a /1 as a suffix to the accession value similarly to the left fragment reads. The right fragment reads should all have /2 as the accession suffix. Then, run Trinity using the --left and --right parameters as if all the data were paired.

*When you have a mix of SE and PE files: for the --left flag, first list the R1 PE files followed by the SE files, then for the --right flag list the R2 PE files*

```
nano trinity_Sid_PE-150_SE-150_control-pop.sh
```

```
#!/bin/bash -l

#SBATCH -o trinity_Sid_PE-150_SE-150_control-pop.txt
#SBATCH -c 32
#SBATCH -p himem
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=trinity_Sid_PE-150_SE-150_control-pop

enable_lmod
module load container_env trinity

crun Trinity --seqType fq --max_memory 768G --no_bowtie --output ./trinity_out_PE-150_SE-150_control-pop --left NC_377A_C_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq,NC_378A_La_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq,NC_379A_N_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq,PC_373A_N_Sid_yr1_R1_clippedtrimmed_nofilter_stillpaired.fastq,PC_365A_C_Sid_yr1_R1_clippedtrimmed_nofilter.fastq,PC_366A_La_Sid_yr1_R1_clippedtrimmed_nofilter.fastq,PC_350A_La_Sid_yr1_R1_clippedtrimmed_nofilter.fastq,PC_364A_N_Sid_yr1_R1_clippedtrimmed_nofilter.fastq,PC_368A_C_Sid_yr1_R1_clippedtrimmed_nofilter.fastq --right NC_377A_C_Sid_yr1_R2_clippedtrimmed_nofilter_stillpaired.fastq,NC_378A_La_Sid_yr1_R2_clippedtrimmed_nofilter_stillpaired.fastq,NC_379A_N_Sid_yr1_R2_clippedtrimmed_nofilter_stillpaired.fastq,PC_373A_N_Sid_yr1_R2_clippedtrimmed_nofilter_stillpaired.fastq --CPU 32
```

```
sbatch trinity_Sid_PE-150_SE-150_control-pop.sh
```

### Check assembly quality

**Best practice to rename with contig number**

grep -c ">" Trinity.fasta
834637

cp Trinity.fasta 834637_Trinity_Sid_PE150-SE150_Control.fasta

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Trinity.fasta
```

```
The total number of sequences is 834637
The average sequence length is 447
The total number of bases is 373170459
The minimum sequence length is 177
The maximum sequence length is 32631
The N50 is 472
Median Length = 565
contigs < 150bp = 0
contigs >= 500bp = 200238
contigs >= 1000bp = 41011
contigs >= 2000bp = 7686
```

### Assess the quality of the Trinity.fasta output 

```
enable_lmod
module load container_env trinity
crun TrinityStats.pl Trinity.fasta
```

Output:

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	690120
Total trinity transcripts:	834637
Percent GC: 41.29

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 1506
	Contig N20: 942
	Contig N30: 707
	Contig N40: 568
	Contig N50: 472

	Median contig length: 340
	Average contig: 447.11
	Total assembled bases: 373170459


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 1171
	Contig N20: 793
	Contig N30: 619
	Contig N40: 511
	Contig N50: 433

	Median contig length: 330
	Average contig: 417.67
	Total assembled bases: 288240225


----------------------------------------------------------------------------------------------
## Longest Isoform Control population assembly

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

nano trinity_Control_LongestIsoform.sh

```
#!/bin/bash -l

#SBATCH -o trinity_Control_LongestIsoform.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=Trinity_Control_LongestIsoform

/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_longest_isoform_per_gene_trinity.py -i Trinity.fasta -o Trinity_Control_LongestIsoform.fasta
```

sbatch trinity_Control_LongestIsoform.sh


**Best practice to rename with contig number**

grep -c ">" Trinity_Control_LongestIsoform.fasta
690120

cp Trinity_Control_LongestIsoform.fasta 690120_Trinity_Sid_Control_LongestIsoform.fasta

#### For some reason, double >> now precede Trinity 'gene' name/ID rather than the single > as before

head -3 Trinity_Control_LongestIsoform.fasta

```
>>TRINITY_DN69561_c0_g1_i1 len=214 path=[0:0-213]
CCAGCTTGGGAATCGGCGGGCGAGATGACAGAGACTGAATTGCGGGGAAGAATAGGTGAATTTTTTATTGTGTTTTCTACTTGTGAGTTCCTCTTTCTTTGCTACATTAGTGCAAGGCTTTTGACTCATAAGCTACAGAAATTAAACAAAGTGTGTGTGATTTAGCAGATATGCTGAATATGCGTGTATTTAGGTGAAGCAATTACAGGCGAAT
>>TRINITY_DN69561_c1_g1_i1 len=279 path=[0:0-278]
```

grep ">>" Trinity_Control_LongestIsoform.fasta | head

Remove first character ('>')
grep ">" Trinity_Control_LongestIsoform.fasta | head | cut -c 2-

cat Trinity_Control_LongestIsoform.fasta | head | cut -c 2- 

```
cat Trinity_Control_LongestIsoform.fasta | cut -c 2- > Trinity_Control_LongestIsoform_edited.fasta
```

head -3 Trinity_Control_LongestIsoform_edited.fasta

grep -c ">" Trinity_Control_LongestIsoform.fasta
690120

grep -c "TRINITY_" Trinity_Control_LongestIsoform.fasta
690120

wc -l Trinity_Control_LongestIsoform.fasta
1380240

*Edited fasta file with single > preceding Trinity sequence/ 'gene' name/ID*

grep -c ">" Trinity_Control_LongestIsoform_edited.fasta
690120

grep -c "TRINITY_" Trinity_Control_LongestIsoform_edited.fasta
690120

wc -l Trinity_Control_LongestIsoform_edited.fasta
1380240

**Best practice to rename with contig number**
cp Trinity_Control_LongestIsoform_edited.fasta 690120_Trinity_Sid_Control_LongestIsoform_edited.fasta

enable_lmod
module load container_env trinity
crun TrinityStats.pl 690120_Trinity_Sid_Control_LongestIsoform_edited.fasta

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	690120
Total trinity transcripts:	690120
Percent GC: 41.11

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 1172
	Contig N20: 793
	Contig N30: 619
	Contig N40: 510
	Contig N50: 432

	Median contig length: 329
	Average contig: 416.67
	Total assembled bases: 287550105


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 1172
	Contig N20: 793
	Contig N30: 619
	Contig N40: 510
	Contig N50: 432

	Median contig length: 329
	Average contig: 416.67
	Total assembled bases: 287550105


### Filter by a length threshold (> 500 bp)

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 500lnThresh 690120_Trinity_Sid_Control_LongestIsoform.fasta
```

Output:
- File:  690120_Trinity_Sid_Control_LongestIsoform500lnThresh.fasta
- Number of total seqs for 690120_Trinity_Sid_Control_LongestIsoform.fasta: 690120
- Number of seqs over 500 for 690120_Trinity_Sid_Control_LongestIsoform.fasta: 146543

**Best practice to rename with contig number**

grep -c ">" 690120_Trinity_Sid_Control_LongestIsoform500lnThresh.fasta
146542

cp 690120_Trinity_Sid_Control_LongestIsoform500lnThresh.fasta 146542_Trinity_Control_LongestIsoform500lnThresh.fasta

### Check the >500 bp contig longest isoform assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 146542_Trinity_Control_LongestIsoform500lnThresh.fasta
```

The total number of sequences is 146542
The average sequence length is 811
The total number of bases is 118941319
The minimum sequence length is 500
The maximum sequence length is 32631
The N50 is 779
Median Length = 729
contigs < 150bp = 0
contigs >= 500bp = 146542
contigs >= 1000bp = 24602
contigs >= 2000bp = 3560


----------------------------------------------------------------------------------------------
## Blast Control population Trinity assembly against Silva rRNA databases (LSU and SSU) to identify rRNA contamination

### BLAST: Basic Local Alignment Search Tool
BLAST finds regions of local similarity between sequences. The program compares nucleotide (or protein) sequences to sequence databases and calculates the statistical significance of matches. Given one or more query sequences (usually in FASTA format), BLAST looks for matching sequence regions between them and a subject set.

A sufficiently close match between subsequences (denoted by arrows in the figure above, though matches are usually longer than illustrated here) is called a high-scoring pair (HSP), while a query sequence is said to hit a target sequence if they share one or more HSPs. Sometimes, however, the term “hit” is used loosely, without differentiating between the two. Each HSP is associated with a “bitscore” that is based on the similarity of the subsequences as determined by a particular set of rules. Because in larger subject sets some good matches are likely to be found by chance, each HSP is also associated with an “E value,” representing the expected number of matches one might find by chance in a subject set of that size with that score or better. For example, an E value of 0.05 means that we can expect a match by chance in 1 in 20 similar searches, whereas an E value of 2.0 means we can expect 2 matches by chance for each similar search.

*Citation:*
Altschul et al. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
[https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=References](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=References)


### SILVA rRNA database 
A comprehensive on-line resource for quality checked and aligned ribosomal RNA sequence data.

SILVA provides comprehensive, quality checked and regularly updated datasets of aligned small (16S/18S, SSU) and large subunit (23S/28S, LSU) ribosomal RNA (rRNA) sequences for all three domains of life (Bacteria, Archaea and Eukarya). 

*Citation:*
Quast et al. (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596. 

Check for latest SILVA database:
[https://www.arb-silva.de/no_cache/download/archive/](https://www.arb-silva.de/no_cache/download/archive/)

Here, I'm using SILVA 138.1 LSU and SSU (August 2020 release).
Access date:  November 5, 2020. 

##### Download directly to cluster

> /cm/shared/courses/dbarshis/barshislab/blastdb/

First file (used wget):

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz

or

curl -O https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz

Second file: 

curl -O https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSUParc_tax_silva_trunc.fasta.gz


##### Uncompress files

gzip -d SILVA_138.1_LSUParc_tax_silva_trunc.fasta.gz

gunzip SILVA_138.1_SSUParc_tax_silva_trunc.fasta.gz


#### Format the LSU and SSU databases for Blast

It turns out that from a computational perspective, simple FASTA files are not easily searched. Thus BLAST+ provides a tool called makeblastdb that converts a subject FASTA file into an indexed and quickly searchable (but not human-readable) version of the same information, stored in a set of similarly named files (often at least three ending in .nin, .nsq, and .nhr for nucleotide sequences). This set of files represents the “database,” and the database name is the shared file name prefix of these files.

Running makeblastdb on a FASTA file: 
```
makeblastdb -in <fasta file> -out <database name> -dbtype <type> -title <title>
```
where <type> is one of prot or nucl, and <title> is a human-readable title (enclosed in quotes if necessary).

blast/2.9
```
enable_lmod
module load container_env blast
makeblastdb -in SILVA_138.1_LSUParc_tax_silva_trunc.fasta -dbtype nucl -out LSUBlastdb -hash_index
```

```
enable_lmod
module load container_env blast
makeblastdb -in SILVA_138.1_SSUParc_tax_silva_trunc.fasta -dbtype nucl -out SSUBlastdb -hash_index
```


----------------------------------------------------------------------------------------------
## Blast SILVA rRNA database against Trinity assembly

File: Trinity_Control_LongestIsoform.fasta
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

## Blast the Trinity assembly to the LSU SILVA database

**Common parameters:**
```
    -query <fasta file>
        The name (or path) of the FASTA-formatted file to search for as query sequences.
    -db <database name>
        The name of the database to search against (as opposed to using -subject).
    -subject <fasta file>
        The name (or path) of the FASTA-formatted file to search in as subject sequences.
    -evalue <real number>
        Only HSPs with E values smaller than this should be reported. For example: -evalue 0.001 or -evalue 1e-6.
    -outfmt <integer>
        How to format the output. The default, 0, provides a human-readable (but not programmatically parseable) text file. The values 6 and 7 produce tab-separated rows and columns in a text file, with 7 providing explanatory comment lines. Similarly, a value of 10 produces comma-separated output; 11 produces a format that can later be quickly turned into any other with another program called blast_formatter. Options 6, 7, and 10 can be highly configured in terms of what columns are shown.
    -num_threads <integer>
        Use <integer> CPU cores on a multicore system, if they are available.
    -max_target_seqs <integer>
        When the output format is 6, 7, or 10 for each query sequence, only report HSPs for the best <integer> different subject sequences.
    -max_hsps <integer>
        For each query/target pair, only report the best <integer> HSPs.
    -out <output file>
        Write the output to <output file> as opposed to the default of standard output.
```
When using the -db option, the BLAST tools will search for the database files in three locations: (1) the present working directory, (2) your home directory, and (3) the paths specified in the $BLASTDB environment variable.

```
nano LSU_Blast_out.sh
```

```
#!/bin/bash -l

#SBATCH -o LSU_Blast_out.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=LSU_Blast_out

enable_lmod
module load container_env blast
blastn -query 146542_Trinity_Control_LongestIsoform500lnThresh.fasta -db /cm/shared/courses/dbarshis/barshislab/blastdb/LSUBlastdb -out LSU_blastn.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
```

```
sbatch LSU_Blast_out.sh
```

## Blast the Trinity assembly to the SSU SILVA database

```
nano SSU_Blast_out.sh
```

```
#!/bin/bash -l

#SBATCH -o SSU_Blast_out.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=SSU_Blast_out

enable_lmod
module load container_env blast
blastn -query 146542_Trinity_Control_LongestIsoform500lnThresh.fasta -db /cm/shared/courses/dbarshis/barshislab/blastdb/SSUBlastdb -out SSU_blastn.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
```

```
sbatch SSU_Blast_out.sh
```

--------------------------------------------------------------------------------------------
## Parse (separate) the blast output 
Remove reads matching to rRNA.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

module spider biopython/1

biopython version 1.75
python/3.6
```
enable_lmod
module load container_env python
module load biopython/1
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py LSUblastn_parsed.txt blastn LSU_blastn.outfmt5
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py SSUblastn_parsed.txt blastn SSU_blastn.outfmt5
```

### Identify "good hits" to SILVA rRNA database
"good hits" are defined here as matching at least 78% of the read over at least 100 bp of the read

NB: original script ReParseBlastbycutoffs_advbioinf.py is for Python 2.
ReParseBlastbycutoffs_advbioinf_b.py edited for Python 3. 

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py _78per100bp.txt lengthidentity 0.78 100 *_parsed.txt
```

```
LSUblastn_parsed.txt	0.78	100	Number of Good Hits:	336
LSUblastn_parsed.txt	0.78	100	Number of unique matches:	285
SSUblastn_parsed.txt	0.78	100	Number of Good Hits:	268
SSUblastn_parsed.txt	0.78	100	Number of unique matches:	191
```

### Remove good LSU and SSU SILVA database hits from the Trinity assembly
Use the getblasthits.py script.

Alternate script (Hannah used):  getseqsfromfasta_advbioinf.py

#### Remove header line
```
tail -n +2 LSUblastn_parsed__78per100bp.txt > LSUblastn_parsed_78per100bp_nohead.txt
tail -n +2 SSUblastn_parsed__78per100bp.txt > SSUblastn_parsed_78per100bp_nohead.txt
```

```
python2 /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py LSUblastn_parsed_78per100bp_nohead.txt 146542_Trinity_Control_LongestIsoform500lnThresh.fasta Trinity_Control_Minus_LSU.fasta
```

```
grep -c ">" Trinity_Control_Minus_LSU.fasta
```
146206

```
python2 /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py SSUblastn_parsed_78per100bp_nohead.txt Trinity_Control_Minus_LSU.fasta Trinity_Control_Minus_LSUSSU.fasta
```

```
grep -c ">" Trinity_Control_Minus_LSUSSU.fasta
```
146151


### Rename the new assembly without rRNA

cp Trinity_Control_Minus_LSUSSU.fasta 146151_Trinity_Control_LongIso-500_MinusSilva.fasta

## Check the new assembly without rRNA

```
python2 /cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 146151_Trinity_Control_LongIso-500_MinusSilva.fasta
```

```
The total number of sequences is 146151
The average sequence length is 811
The total number of bases is 118637227
The minimum sequence length is 500
The maximum sequence length is 32631
The N50 is 779
Median Length = 1118
contigs < 150bp = 0
contigs >= 500bp = 146151
contigs >= 1000bp = 24556
contigs >= 2000bp = 3549
```


#### Edit >> to single > at the beginning of each Trinity 'gene' name ()
For some reason double >> appeared again?
```
head -3 146151_Trinity_Control_LongIso-500_MinusSilva.fasta
```

*Need to remove first character (the extra '>')*
grep ">" 146151_Trinity_Control_LongIso-500_MinusSilva.fasta | head
grep ">" 146151_Trinity_Control_LongIso-500_MinusSilva.fasta | head | cut -c 2-
cat 146151_Trinity_Control_LongIso-500_MinusSilva.fasta | head | cut -c 2- 

```
cat 146151_Trinity_Control_LongIso-500_MinusSilva.fasta | cut -c 2- > 146151_Trinity_Control_LongIso-500_MinusSilva_edited.fasta
```

**Best practice to rename with contig number**

grep -c ">" 146151_Trinity_Control_LongIso-500_MinusSilva_edited.fasta
146151

cp 146151_Trinity_Control_LongIso-500_MinusSilva_edited.fasta 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta


----------------------------------------------------------------------------------------------
## Separation of Coral Host and Symbiont contigs

Use the Control population, longest isoform, >500 bp contigs minus SILVA assembly.
Trinity assembled reference with 
- longest isoforms
- rRNA removed 
- only contigs greater than 500 bp

Same as file:  146151_Trinity_Control_LongIso-500_Minus_Silva.fasta

#### Create Blast databases from Hannah's clean coral (CC), dirty coral (DC), clean symbiont (CS) and dirty symbiont (DS) references

Location of Hannah's databases:
> /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/

```
enable_lmod
module load container_env blast
makeblastdb -in /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CC_all_db.fasta -dbtype nucl -out CC_all_Blastdb -hash_index
```

```
Building a new DB, current time: 01/14/2021 16:11:35
New DB name:   /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/CC_all_Blastdb
New DB title:  /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CC_all_db.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Ignoring sequence 'lcl|4841' as it has no sequence data
Ignoring sequence 'lcl|475475' as it has no sequence data
Ignoring sequence 'lcl|475476' as it has no sequence data
Adding sequences from FASTA; added 847958 sequences in 107.903 seconds.
```

Created 8 files:

CC_all_Blastdb.nhd
CC_all_Blastdb.nhi
CC_all_Blastdb.nhr
CC_all_Blastdb.nin
CC_all_Blastdb.nog
CC_all_Blastdb.nsd
CC_all_Blastdb.nsi
CC_all_Blastdb.nsq


Re-make the rest of the databases:
```
makeblastdb -in /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CS_all_db.fasta -dbtype nucl -out CS_all_Blastdb -hash_index
```

Building a new DB, current time: 01/14/2021 16:19:29
New DB name:   /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/CS_all_Blastdb
New DB title:  /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CS_all_db.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
FASTA-Reader: Ignoring invalid residues at position(s): On line 3871821: 57-58
*then a gazillion more lines: "FASTA-Reader: Ignoring invalid residues at position(s): On line ..."*
Adding sequences from FASTA; added 403809 sequences in 33.7596 seconds.

^ Created 8 CS files.

```
makeblastdb -in /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/DC_all_db.fasta -dbtype nucl -out DC_all_Blastdb -hash_index
```

Building a new DB, current time: 01/14/2021 16:23:12
New DB name:   /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/DC_all_Blastdb
New DB title:  /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/DC_all_db.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
FASTA-Reader: Ignoring invalid residues at position(s): On line 66992: 429-432, 434-438
FASTA-Reader: Ignoring invalid residues at position(s): On line 9045440: 1119-1122, 1124-1125, 1127-1128, 1130-1132, 1135-1136, 1141-1143, 1145, 1147-1156, 1159, 1161-1162, 1164-1165, 1170, 1172-1192
Ignoring sequence 'lcl|2235324' as it has no sequence data
Ignoring sequence 'lcl|2705958' as it has no sequence data
Ignoring sequence 'lcl|2705959' as it has no sequence data
Adding sequences from FASTA; added 3078441 sequences in 266.807 seconds.

^ Created 17 DC files.

```
makeblastdb -in /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/DS_all_db.fasta -dbtype nucl -out DS_all_Blastdb -hash_index
```

Building a new DB, current time: 01/14/2021 16:19:29
New DB name:   /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/CS_all_Blastdb
New DB title:  /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CS_all_db.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
FASTA-Reader: Ignoring invalid residues at position(s): On line 3871821: 57-58
*then a gazillion more lines: "FASTA-Reader: Ignoring invalid residues at position(s): On line ..."*
Adding sequences from FASTA; added 747502 sequences in 61.2529 seconds.

^ Created 8 DS files. 



#### Blast the file to the clean coral (CC), dirty coral (DC), clean symbiont (CS) and dirty symbiont (DS) databases (Hannah's)

```
nano blastTrinitytoAnnotate.sh
```

```
#!/bin/bash -l

#SBATCH -o blastTrinitytoAnnotate.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=blastTrinitytoAnnotate

enable_lmod
module load container_env blast

blastn -query Trinity_Control_LongIso-500_Minus_Silva.fasta -db /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CC_all_nucldb -out TrinityMinusSilva_CC.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
blastn -query Trinity_Control_LongIso-500_Minus_Silva.fasta -db /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/DC_all_nucldb -out TrinityMinusSilva_DC.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
blastn -query Trinity_Control_LongIso-500_Minus_Silva.fasta -db /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CS_all_nucldb -out TrinityMinusSilva_CS.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
blastn -query Trinity_Control_LongIso-500_Minus_Silva.fasta -db /cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/DS_all_nucldb -out TrinityMinusSilva_DS.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
```

```
sbatch blastTrinitytoAnnotate.sh
```

#### Check Blast outputs
The Blast outputs are the *.outfmt5 files.

grep "<Iteration_query-len>" TrinityMinusSilva_CC.outfmt5 | wc -l
146151

grep "<Iteration_query-len>" TrinityMinusSilva_CS.outfmt5 | wc -l
146151

grep "<Iteration_query-len>" TrinityMinusSilva_DC.outfmt5 | wc -l
146151

grep "<Iteration_query-len>" TrinityMinusSilva_DS.outfmt5 | wc -l
146151

```
head TrinityMinusSilva_CC.outfmt5
```

```
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.9.0+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), &quot;A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db>/cm/shared/courses/dbarshis/RedSea/blastdbs/New_blast_dbs/CC_all_nucldb</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>&gt;TRINITY_DN521319_c0_g1_i2 len=804 path=[0:0-803]</BlastOutput_query-def>
  <BlastOutput_query-len>804</BlastOutput_query-len>
```

A few fields appear to specify the contig length. Figure out a way to isolate just the length:
```
grep "<Iteration_query-len>" TrinityMinusSilva_CC.outfmt5 | head | cut -d ">" -f 2
grep "<Iteration_query-len>" TrinityMinusSilva_CC.outfmt5 | head | cut -d ">" -f 2 | cut -d "<" -f 1
grep "<Iteration_query-len>" TrinityMinusSilva_CC.outfmt5 | cut -d ">" -f 2 | cut -d "<" -f 1 | sort -n | head

grep "<Iteration_query-len>" TrinityMinusSilva_CS.outfmt5 | head | cut -d ">" -f 2 | cut -d "<" -f 1
grep "<Iteration_query-len>" TrinityMinusSilva_CS.outfmt5 | cut -d ">" -f 2 | cut -d "<" -f 1 | sort -n | head

grep "<Iteration_query-len>" TrinityMinusSilva_DC.outfmt5 | cut -d ">" -f 2 | cut -d "<" -f 1 | sort -n | head

grep "<Iteration_query-len>" TrinityMinusSilva_DS.outfmt5 | cut -d ">" -f 2 | cut -d "<" -f 1 | sort -n | head
```

**^ Confirmed that all contig lengths are a minimum of 500 bp.**

#### Check .fasta files
Check output for the post-blast .fasta files versus pre-blast .fasta files

```
bash
for i in *.fasta ; do echo $i; /cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py $i ; done
```

*The TRINITY 'gene' names are on different lines than their respective sequences*

*In the post-blast files, the Trinity 'gene' names are on the SAME line as the sequence and other miscellaneous Blast information*

head -1 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta
```
>TRINITY_DN521319_c0_g1_i2 len=804 path=[0:0-803]
```

head -1 11704_Trinity-Silva_CC_parsed_70per100bpmatch.fasta
```
>TRINITY_DN22892_c0_g1_i2 len=638 path=[0:0-445 2:446-637]	638	gnl|BL_ORD_ID|411914 s23Contig8188_9 [1660 - 2211] 	552	346	288	633	551	206	AGCCATATGTCTGAAGAACAGTCTCCCCCCGGATATCGCATTTCATGAGCAAGTGCAGGACCATCAGGCTCCTCTCCAAAGTCCACGCAGAAATCAGGATTCAACTTCAGACCACCGTTCACGGAATCTACGTCCACTTGGAGCAACATGCTGCCTTTCTTGGACATATTTGGGTAGAACTGATTGTCCCACACGCTGTACAAAGATGTGGTAACGTACAGTCGTTTTCCATCCAAACTAAGCTGAATCATTTGCGGTCCTCCTTCCACGCGCTTCCCCTTCATATAACACGGCTCAGGCTGTGCACTCAGTTCAGAGTCCTCGATCACCTTCACTGGGCCATCAC	AGCCATATATCAGAGGAGCAGTCTCCTCCCGGATAACGCACTTCATGAGCGAGTGCTGGTCCGTCAGGTTCTTCTCCAAAGTCCACACAGAAATCAGGATTCAACTTGAGCCCACCATTCACTGAATCAACGTCAACTTGGAGCAACATGCTGCCTTTCTTGGAGAGATTAGGGTAGAACTGATTGTCCCACACACTGAACAACGATGTTGTGACGTAAAGCCGTTTGCCATCCAGACTGAGCTGAATCATCTGCGGTCCTCCTTCCACTCGCTTCCCCTTCACGTAACACGGCTCAGGCTGTCCACTCAGTTCAGAGTCCTCTGTTACCTTGACTGGGCCATCAC	220.0	5.49555e-111	304	0.8786127167630058	0
```


----------------------------------------------------------------------------------------------
## Parse the blast outputs to each of the four databases using the same parse_blastnortblastx_advbioinf.py script

```
enable_lmod
module load container_env python
module load biopython/1

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py TrinityMinusSilva_CC_parsed.txt blastn TrinityMinusSilva_CC.outfmt5

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py TrinityMinusSilva_CS_parsed.txt blastn TrinityMinusSilva_CS.outfmt5

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py TrinityMinusSilva_DC_parsed.txt blastn TrinityMinusSilva_DC.outfmt5

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py TrinityMinusSilva_DS_parsed.txt blastn TrinityMinusSilva_DS.outfmt5
```

----------------------------------------------------------------------------------------------
## Re-Parse the blast outputs 

Use the ReParseBlastbycutoffs.py script with various cut-offs for 
- the percentage of the read that has to match (first number input after 'length identity') 
- the length of the read it has to match (second number input after 'length identity'). 
These different options will be utilized in the next step.

NB: original script ReParseBlastbycutoffs_advbioinf.py is for Python 2.
ReParseBlastbycutoffs_advbioinf_b.py edited for Python 3. 

### Identify "good hits" to Hannah's assembly
"good hits" are defined here as matching at least __% of the read over at least 100 bp of the read

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 60per100bpmatch.txt lengthidentity 0.60 100 TrinityMinusSilva_*_parsed.txt
```

```
TrinityMinusSilva_CC_parsed.txt	0.60	100	Number of Good Hits:	11704
TrinityMinusSilva_CC_parsed.txt	0.60	100	Number of unique matches:	7488
TrinityMinusSilva_CS_parsed.txt	0.60	100	Number of Good Hits:	3006
TrinityMinusSilva_CS_parsed.txt	0.60	100	Number of unique matches:	2690
TrinityMinusSilva_DC_parsed.txt	0.60	100	Number of Good Hits:	26354
TrinityMinusSilva_DC_parsed.txt	0.60	100	Number of unique matches:	21580
TrinityMinusSilva_DS_parsed.txt	0.60	100	Number of Good Hits:	8750
TrinityMinusSilva_DS_parsed.txt	0.60	100	Number of unique matches:	7537
```

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 70per100bpmatch.txt lengthidentity 0.70 100 TrinityMinusSilva_*_parsed.txt
```

```
TrinityMinusSilva_CC_parsed.txt	0.70	100	Number of Good Hits:	11704
TrinityMinusSilva_CC_parsed.txt	0.70	100	Number of unique matches:	7488
TrinityMinusSilva_CS_parsed.txt	0.70	100	Number of Good Hits:	3006
TrinityMinusSilva_CS_parsed.txt	0.70	100	Number of unique matches:	2690
TrinityMinusSilva_DC_parsed.txt	0.70	100	Number of Good Hits:	26352
TrinityMinusSilva_DC_parsed.txt	0.70	100	Number of unique matches:	21578
TrinityMinusSilva_DS_parsed.txt	0.70	100	Number of Good Hits:	8750
TrinityMinusSilva_DS_parsed.txt	0.70	100	Number of unique matches:	7537
```

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 75per100bpmatch.txt lengthidentity 0.75 100 TrinityMinusSilva_*_parsed.txt
```

```
TrinityMinusSilva_CC_parsed.txt	0.75	100	Number of Good Hits:	11350
TrinityMinusSilva_CC_parsed.txt	0.75	100	Number of unique matches:	7321
TrinityMinusSilva_CS_parsed.txt	0.75	100	Number of Good Hits:	2894
TrinityMinusSilva_CS_parsed.txt	0.75	100	Number of unique matches:	2584
TrinityMinusSilva_DC_parsed.txt	0.75	100	Number of Good Hits:	25825
TrinityMinusSilva_DC_parsed.txt	0.75	100	Number of unique matches:	21204
TrinityMinusSilva_DS_parsed.txt	0.75	100	Number of Good Hits:	8731
TrinityMinusSilva_DS_parsed.txt	0.75	100	Number of unique matches:	7520
```

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 80per100bpmatch.txt lengthidentity 0.80 100 TrinityMinusSilva_*_parsed.txt
```

```
TrinityMinusSilva_CC_parsed.txt	0.80	100	Number of Good Hits:	9359
TrinityMinusSilva_CC_parsed.txt	0.80	100	Number of unique matches:	6322
TrinityMinusSilva_CS_parsed.txt	0.80	100	Number of Good Hits:	2301
TrinityMinusSilva_CS_parsed.txt	0.80	100	Number of unique matches:	2059
TrinityMinusSilva_DC_parsed.txt	0.80	100	Number of Good Hits:	22456
TrinityMinusSilva_DC_parsed.txt	0.80	100	Number of unique matches:	18701
TrinityMinusSilva_DS_parsed.txt	0.80	100	Number of Good Hits:	8546
TrinityMinusSilva_DS_parsed.txt	0.80	100	Number of unique matches:	7366
```

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 90per100bpmatch.txt lengthidentity 0.90 100 TrinityMinusSilva_*_parsed.txt
```

```
TrinityMinusSilva_CC_parsed.txt	0.90	100	Number of Good Hits:	1600
TrinityMinusSilva_CC_parsed.txt	0.90	100	Number of unique matches:	1305
TrinityMinusSilva_CS_parsed.txt	0.90	100	Number of Good Hits:	131
TrinityMinusSilva_CS_parsed.txt	0.90	100	Number of unique matches:	122
TrinityMinusSilva_DC_parsed.txt	0.90	100	Number of Good Hits:	9420
TrinityMinusSilva_DC_parsed.txt	0.90	100	Number of unique matches:	8365
TrinityMinusSilva_DS_parsed.txt	0.90	100	Number of Good Hits:	7423
TrinityMinusSilva_DS_parsed.txt	0.90	100	Number of unique matches:	6535
```


----------------------------------------------------------------------------------------------
## Extract the 70 bp per 100 match parsed hits from the Trinity assembly

#### Remove header line
```
tail -n +2 TrinityMinusSilva_CC_parsed_70per100bpmatch.txt > TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt
tail -n +2 TrinityMinusSilva_CS_parsed_70per100bpmatch.txt > TrinityMinusSilva_CS_parsed_70per100bpmatch_nohead.txt
tail -n +2 TrinityMinusSilva_DC_parsed_70per100bpmatch.txt > TrinityMinusSilva_DC_parsed_70per100bpmatch_nohead.txt
tail -n +2 TrinityMinusSilva_DS_parsed_70per100bpmatch.txt > TrinityMinusSilva_DS_parsed_70per100bpmatch_nohead.txt
```

### Make HITS list with the contig name (get Blast hits) and make into fasta files
Use getblasthits.py script

head -1 TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt
```
>TRINITY_DN22892_c0_g1_i2 len=638 path=[0:0-445 2:446-637]	638	gnl|BL_ORD_ID|411914 s23Contig8188_9 [1660 - 2211] 	552	346	288	633	551	206	AGCCATATGTCTGAAGAACAGTCTCCCCCCGGATATCGCATTTCATGAGCAAGTGCAGGACCATCAGGCTCCTCTCCAAAGTCCACGCAGAAATCAGGATTCAACTTCAGACCACCGTTCACGGAATCTACGTCCACTTGGAGCAACATGCTGCCTTTCTTGGACATATTTGGGTAGAACTGATTGTCCCACACGCTGTACAAAGATGTGGTAACGTACAGTCGTTTTCCATCCAAACTAAGCTGAATCATTTGCGGTCCTCCTTCCACGCGCTTCCCCTTCATATAACACGGCTCAGGCTGTGCACTCAGTTCAGAGTCCTCGATCACCTTCACTGGGCCATCAC	AGCCATATATCAGAGGAGCAGTCTCCTCCCGGATAACGCACTTCATGAGCGAGTGCTGGTCCGTCAGGTTCTTCTCCAAAGTCCACACAGAAATCAGGATTCAACTTGAGCCCACCATTCACTGAATCAACGTCAACTTGGAGCAACATGCTGCCTTTCTTGGAGAGATTAGGGTAGAACTGATTGTCCCACACACTGAACAACGATGTTGTGACGTAAAGCCGTTTGCCATCCAGACTGAGCTGAATCATCTGCGGTCCTCCTTCCACTCGCTTCCCCTTCACGTAACACGGCTCAGGCTGTCCACTCAGTTCAGAGTCCTCTGTTACCTTGACTGGGCCATCAC	220.0	5.49555e-111	304	0.8786127167630058	0
```
*May have single or multiple paths - variable!*

```
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt | head
```
TRINITY_DN22892_c0_g1_i2
TRINITY_DN64979_c0_g1_i6
TRINITY_DN484077_c0_g1_i1

```
grep -o -E "(TRINITY_\w+\s\w+.\w+)" TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt | head
```
TRINITY_DN22892_c0_g1_i2 len=638
TRINITY_DN64979_c0_g1_i6 len=578
TRINITY_DN484077_c0_g1_i1 len=531

***Linux server grep likely uses Perl Compatible Regular Expressions (PCRE) regex engine*** 
-P flag for Perl
```
enable_lmod
module load container_env perl
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt | head
```

**I DID IT!**
```
TRINITY_DN22892_c0_g1_i2 len=638 path=[0:0-445 2:446-637]
TRINITY_DN64979_c0_g1_i6 len=578 path=[2:0-110 4:111-208 7:209-242 8:243-256 9:257-369 11:370-577]
TRINITY_DN484077_c0_g1_i1 len=531 path=[0:0-530]
TRINITY_DN63583_c0_g1_i3 len=1510 path=[0:0-425 1:426-653 2:654-702 4:703-788 5:789-1509]
TRINITY_DN91634_c0_g1_i1 len=551 path=[1:0-95 2:96-550]
TRINITY_DN37872_c0_g1_i2 len=828 path=[0:0-583 1:584-678 3:679-827]
comp90172_c0_seq1 len=671 path=[1:0-238 240:239-670]
TRINITY_DN36661_c0_g1_i1 len=1412 path=[1:0-1225 2:1226-1411]
TRINITY_DN551315_c0_g1_i1 len=706 path=[0:0-705]
TRINITY_DN35809_c6_g1_i1 len=567 path=[0:0-129 2:130-566]
```

### Create HITS file
- Clean Coral (CC)
```
enable_lmod
module load container_env perl
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" TrinityMinusSilva_CC_parsed_70per100bpmatch_nohead.txt > TrinityMinusSilva_CC_parsed_70per100bpmatch_HITS.txt
```

```
head -3 TrinityMinusSilva_CC_parsed_70per100bpmatch_HITS.txt
```
TRINITY_DN22892_c0_g1_i2 len=638 path=[0:0-445 2:446-637]
TRINITY_DN64979_c0_g1_i6 len=578 path=[2:0-110 4:111-208 7:209-242 8:243-256 9:257-369 11:370-577]
TRINITY_DN484077_c0_g1_i1 len=531 path=[0:0-530]

- Clean Symbiont
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" TrinityMinusSilva_CS_parsed_70per100bpmatch_nohead.txt > TrinityMinusSilva_CS_parsed_70per100bpmatch_HITS.txt

- Dirty Coral
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" TrinityMinusSilva_DC_parsed_70per100bpmatch_nohead.txt > TrinityMinusSilva_DC_parsed_70per100bpmatch_HITS.txt

- Dirty Symbiont
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" TrinityMinusSilva_DS_parsed_70per100bpmatch_nohead.txt > TrinityMinusSilva_DS_parsed_70per100bpmatch_HITS.txt


### Get blast hits from each database (CC, CS, DC, DS) and create individual fasta files

##### Clean Coral database
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py TrinityMinusSilva_CC_parsed_70per100bpmatch_HITS.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta TrinityMinusSilva_CC_parsed_70per100bpmatch.fasta
```

```
grep -c ">" TrinityMinusSilva_CC_parsed_70per100bpmatch.fasta
```
134447

Alternate- If Python 3 already loaded:
```
enable_lmod
module load container_env python
python2 /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py TrinityMinusSilva_CC_parsed_70per100bpmatch_HITS.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta TrinityMinusSilva_CC_parsed_70per100bpmatch.fasta
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Trinity_Control_LongIso-500_MinusSilva_CC_70-100match.fasta
```

```
The total number of sequences is 134447
The average sequence length is 776
The total number of bases is 104441196
The minimum sequence length is 499
The maximum sequence length is 15410
The N50 is 748
Median Length = 569
contigs < 150bp = 0
contigs >= 500bp = 133810
contigs >= 1000bp = 20100
contigs >= 2000bp = 2190
```

##### Clean Symbiont database
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py TrinityMinusSilva_CS_parsed_70per100bpmatch_HITS.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta TrinityMinusSilva_CS_parsed_70per100bpmatch.fasta
```

```
grep -c ">" TrinityMinusSilva_CS_parsed_70per100bpmatch.fasta
```
143145

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py TrinityMinusSilva_CS_parsed_70per100bpmatch.fasta
```

```
The total number of sequences is 143145
The average sequence length is 809
The total number of bases is 115895601
The minimum sequence length is 499
The maximum sequence length is 32630
The N50 is 777
Median Length = 569
contigs < 150bp = 0
contigs >= 500bp = 142496
contigs >= 1000bp = 23807
contigs >= 2000bp = 3456
```

##### Dirty Coral database
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py TrinityMinusSilva_DC_parsed_70per100bpmatch_HITS.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta TrinityMinusSilva_DC_parsed_70per100bpmatch.fasta
```

```
grep -c ">" TrinityMinusSilva_DC_parsed_70per100bpmatch.fasta
```
119799


```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py TrinityMinusSilva_DC_parsed_70per100bpmatch.fasta
```

```
The total number of sequences is 119799
The average sequence length is 761
The total number of bases is 91251200
The minimum sequence length is 499
The maximum sequence length is 15410
The N50 is 734
Median Length = 1271
contigs < 150bp = 0
contigs >= 500bp = 119222
contigs >= 1000bp = 16677
contigs >= 2000bp = 1590
```

##### Dirty Symbiont database
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/getblasthits.py TrinityMinusSilva_DS_parsed_70per100bpmatch_HITS.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta TrinityMinusSilva_DS_parsed_70per100bpmatch.fasta
```

grep -c ">" TrinityMinusSilva_DS_parsed_70per100bpmatch.fasta
137401

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py TrinityMinusSilva_DS_parsed_70per100bpmatch.fasta
```

```
The total number of sequences is 137401
The average sequence length is 808
The total number of bases is 111126036
The minimum sequence length is 499
The maximum sequence length is 32630
The N50 is 776
Median Length = 735
contigs < 150bp = 0
contigs >= 500bp = 136781
contigs >= 1000bp = 22837
contigs >= 2000bp = 3248
```


----------------------------------------------------------------------------------------------
## Taxon assignment - part A (Coral)
### Cross-check & clean Coral and Symbiont assemblies

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

Here we are using 70% identity cutoff matches.

## To get "Good" (clean) Coral [GC]

Start with all the hits >70% id and >100bp overlap with dirty coral (DC) (everything that might be a coral).

### 1. Remove any hits that are also >70% id and >100bp overlap with clean symbiont (anything that is very likely a symbiont)
*Goal:  select all of the “good” coral matches.*

grep -c ">" TrinityMinusSilva_DC_parsed_70per100bpmatch.txt
26352

#### Compare contig names in DC versus CS files 

```
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_DC_parsed_70per100bpmatch.txt > DC_contig.txt
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_CS_parsed_70per100bpmatch.txt > CS_contig.txt
```

Check if any overlapping contigs?
*There are 2903 same TRINITY names (but not necessarily same sequences / blast hits) in CS that are also in DC*

#### Search for & remove duplicate Trinity 'gene' names that exist in CS and DC

*Check for and remove matches between two dataframes*
```
enable_lmod
module load container_env python
module load container_env pandas

python compareDataframes.py TrinityMinusSilva_CS_parsed_70per100bpmatch.txt TrinityMinusSilva_DC_parsed_70per100bpmatch.txt DC_clean.txt CS_clean.txt
```

grep -c ">" DC_clean.txt
23398

grep -c ">" CS_clean.txt
52

**Best practice to rename with contig number**
cp DC_clean.txt 23398_Sid_DCcleaned.txt

### 2. Double check that all of the DC_cleaned are in your original CC (clean coral) matches.

grep -c ">" TrinityMinusSilva_CC_parsed_70per100bpmatch.txt
11704

If you want to cross-check potential duplicates between the files:
```
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_CC_parsed_70per100bpmatch.txt > CC_contig.txt
grep -o -E "(TRINITY_\w+)" DC_clean.txt > DCclean_contig.txt
```

Check if any CC (clean coral) database hits are missing from the DC_cleaned; if CC contigs are missing, add them to the DC_cleaned assembly
```
enable_lmod
module load container_env python
module load container_env pandas

python compareConcatDf.py TrinityMinusSilva_CC_parsed_70per100bpmatch.txt 23398_Sid_DCcleaned.txt DCcleaned_CCchecked.txt
```

grep -c ">" DCcleaned_CCchecked.txt
23842

**Best practice to rename with contig number**
cp DCcleaned_CCchecked.txt 23842_Sid_DCcleaned_CCchecked.txt
cp DCcleaned_CCchecked.txt 23842_Sid_DCcleaned_CCchecked.fasta


### 3. Remove anything that is also in your “good” symbiont reference.

*First you need to make the Good Symbiont reference (next step)*

Double check if any cleaned Good Coral contigs remain in the cleaned Symbiont file 7651_SidSymb_DScleaned_CSchecked.txt
```
enable_lmod
module load container_env python
module load container_env pandas

python compareDataframes.py 7651_SidSymb_DScleaned_CSchecked.txt 23842_Sid_DCcleaned_CCchecked.txt Sid_GoodCoral.txt SidSymb_DScleaned_CSchecked_GCchecked.txt
```

grep -c ">" Sid_GoodCoral.txt
19272

Can cross check this with Good Symbiont output further below. Confirmed below. 
grep -c ">" SidSymb_DScleaned_CSchecked_GCchecked.txt
3081

If you want to confirm potential duplicates between the files:
*4570 Trinity 'gene' names were found in BOTH coral and symbiont files...seems like a lot?*

The overlapping sequences are ones that match things in both the dirtycoral and dirtysym databases but are not found in either the cleancoral or cleansym so we can't definitively say they are host or symbiont. 
It's not a perfect process but it'd be better to toss them out now versus have a potential for a symbiont gene to influence our coral results or vice versa.

```
grep -o -E "(TRINITY_\w+)" 7651_SidSymb_DScleaned_CSchecked.txt > 7651_SidSymb_DScleaned_CSchecked_contig.txt
grep -o -E "(TRINITY_\w+)" 23842_Sid_DCcleaned_CCchecked.txt > 23842_Sid_DCcleaned_CCchecked.txt_contig.txt
```


## Extract good Coral contigs from the Trinity minus SILVA, >500 bp contigs assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/getseqsfromfasta_advbioinf.py 19272_Sid_GoodCoral.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta Sid_GoodCoral_Final.fasta
```

#### Extract sequences with names by creating new dataframe and saving as fasta file
- Query Name
- Query Sequence

```
enable_lmod
module load container_env python
module load container_env pandas
python extractColumns.py Sid_GoodCoral.txt Sid_GoodCoral.fasta
```

grep -c ">" Sid_GoodCoral.fasta
19272

```
enable_lmod
module load container_env perl
grep -P -o "(.\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" Sid_GoodCoral.fasta | head
grep -P -o "(.\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" Sid_GoodCoral.fasta > Sid_GoodCoral_HITS.txt
```

grep -c ">" Sid_GoodCoral_HITS.txt
19272

grep -c ">" Sid_GoodCoral_Final.fasta
19272

Output file:  Sid_GoodCoral_Final.txt

*Could have used the cut command in bash to break up lines/ isolate parts of a line based on delimiters*

**Best practice to rename with contig number**
cp Sid_GoodCoral_Final.txt 19272_Sid_GoodCoral_Final.fasta

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 19272_Sid_GoodCoral_Final.fasta
```

```
The total number of sequences is 19272
The average sequence length is 1125
The total number of bases is 21692610
The minimum sequence length is 499
The maximum sequence length is 32630
The N50 is 1211
Median Length = 582
contigs < 150bp = 0
contigs >= 500bp = 19223
contigs >= 1000bp = 6693
contigs >= 2000bp = 1829
```


----------------------------------------------------------------------------------------------
## Taxon assignment - part B (Symbiont)
### Cross-check & clean Coral and Symbiont assemblies

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/

Here we are using 70% identity cutoff matches.

## To get "Good" (clean) Symbiont [GS]

Start with all the hits >70% id and >100bp overlap with dirty symbiont (DS) (everything that might be a symbiont).

### 1. Remove any hits that are also >70% id and >100bp overlap with clean coral [CC] (anything that is very likely a coral)
Goal:  select all of the “good” symbiont matches.

grep -c ">" TrinityMinusSilva_DS_parsed_70per100bpmatch.txt
8750

#### Compare contig names in DS versus CC files 

```
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_DS_parsed_70per100bpmatch.txt > DS_contig.txt
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_CC_parsed_70per100bpmatch.txt > CC_contig.txt
```

#### Search for & remove duplicate Trinity 'gene' names that exist in DS and CC

grep -c ">" TrinityMinusSilva_CC_parsed_70per100bpmatch.txt
11704

grep -c ">" TrinityMinusSilva_DS_parsed_70per100bpmatch.txt
8750

*Check for and remove matches between two dataframes*
```
enable_lmod
module load container_env python
module load container_env pandas

python compareDataframes.py TrinityMinusSilva_CC_parsed_70per100bpmatch.txt TrinityMinusSilva_DS_parsed_70per100bpmatch.txt DS_clean.txt CC_clean.txt
```

grep -c ">" DS_clean.txt
7301

grep -c ">" CC_clean.txt
10255

**Best practice to rename with contig number**
cp DS_clean.txt 7301_SidSymb_DScleaned.txt


### 2. Double check that all of the DS_cleaned are in your original CS (clean symbiont) matches.

grep -c ">" TrinityMinusSilva_CS_parsed_70per100bpmatch.txt
3006

```
grep -o -E "(TRINITY_\w+)" TrinityMinusSilva_CS_parsed_70per100bpmatch.txt > CS_contig.txt
grep -o -E "(TRINITY_\w+)" DS_clean.txt > DSclean_contig.txt
```

Check if any CS (clean cultured symbiont) database hits are missing from the DS_cleaned; if CS contigs are missing, add them to the DS_cleaned assembly
```
enable_lmod
module load container_env python
module load container_env pandas

python compareConcatDf.py TrinityMinusSilva_CS_parsed_70per100bpmatch.txt DS_clean.txt DScleaned_CSchecked.txt
```

grep -c ">" DScleaned_CSchecked.txt
7651

**Best practice to rename with contig number**
cp DScleaned_CSchecked.txt 7651_SidSymb_DScleaned_CSchecked.txt
cp DScleaned_CSchecked.txt 7651_SidSymb_DScleaned_CSchecked.fasta


### 3. Remove anything that is also in your “good” coral reference.

Double check if any cleaned Good Coral contigs remain in the cleaned Symbiont file 7651_SidSymb_DScleaned_CSchecked.txt
```
enable_lmod
module load container_env python
module load container_env pandas

python compareDataframes.py 23842_Sid_DCcleaned_CCchecked.txt 7651_SidSymb_DScleaned_CSchecked.txt SidSymb_GoodSymb.txt Sid_DCcleaned_CCchecked_GSchecked.txt
```

grep -c ">" SidSymb_GoodSymb.txt
3081

This is consistent with above Good Coral output. 
grep -c ">" Sid_DCcleaned_CCchecked_GSchecked.txt
19272

Triple check (check how many duplicate Trinity 'genes' are present in each of the files)
*4570 Trinity 'gene' names were found in BOTH coral and symbiont files...seems like a lot?*

Consistent with what found above when filtering Good Coral.
```
grep -o -E "(TRINITY_\w+)" 23842_Sid_DCcleaned_CCchecked.txt > Sid_DCcleaned_CCchecked_contig.txt
grep -o -E "(TRINITY_\w+)" 7651_SidSymb_DScleaned_CSchecked.txt > SidSymb_DScleaned_CSchecked_contig.txt
```

As for the symbionts, since we'll be including the genomes anyway it shouldn't be too big of an issue to have that few de novo assembled contigs. 

Regardless we should be a little cautious with analyzing the symbiont data anyway given the low coverage and potentially mixed genera.

(NB: Later we decide not to use symbiont de novo assembly)


## Extract good Symbiont contigs from the Trinity minus SILVA, >500 bp contigs assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/getseqsfromfasta_advbioinf.py 3081_Sid_GoodSymb.txt 146151_Trinity_Control_LongIso-500_Minus_Silva.fasta Sid_GoodSymbiont_Final.fasta
```


#### Extract sequences with names by creating new dataframe and saving as fasta file
- Query Name
- Query Sequence

```
enable_lmod
module load container_env python
module load container_env pandas
python extractColumns.py SidSymb_GoodSymb.txt SidSymb_GoodSymb.fasta
```

grep -c ">" SidSymb_GoodSymb.fasta
3081

```
enable_lmod
module load container_env perl
grep -P -o "(.\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" SidSymb_GoodSymb.fasta | head
grep -P -o "(.\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" SidSymb_GoodSymb.fasta > Sid_GoodSymb_HITS.txt
```

grep -c ">" Sid_GoodSymb_HITS.txt
3081

grep -c ">" Sid_GoodSymb_Final.txt
3081

*Could have used the cut command in bash to break up lines/ isolate parts of a line based on delimiters*

**Best practice to rename with contig number**
cp Sid_GoodSymb_Final.txt 3081_Sid_GoodSymb_Final.fasta

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 3081_Sid_GoodSymb_Final.fasta
```

```
The total number of sequences is 3081
The average sequence length is 842
The total number of bases is 2594909
The minimum sequence length is 499
The maximum sequence length is 6438
The N50 is 841
Median Length = 570
contigs < 150bp = 0
contigs >= 500bp = 3063
contigs >= 1000bp = 657
contigs >= 2000bp = 81
```


----------------------------------------------------------------------------------------------
## Rename Good Coral and Good Symbiont references

### Coral Host assembly:
19272_Sid_GoodCoral_Final.fasta


### Symbiont assembly:
3081_Sid_GoodSymb_Final.fasta

```
head -2 3081_Sid_GoodSymb_Final.fasta
```

```
>TRINITY_DN454986_c0_g1_i1 len=621 path=[0:0-620]
GGTTGGTGTCAAGCATCATGCCCGACGTTTGCAAGCCCTAGAGGGAGTTTCAACAAGTGACATTATTTCTTCTTTGCCCCGCGAAGAGCTGGTCGCGAAGATCCGCAGCGCTTTGCCATACCTTGCCCAACTGGCAAAGGGCGGCAGTCGCCCCACTGCAGCACTTGAGAAGTCTATGATTGACACGACCACAGCGATGCCACAGGCCAAAGCCACCGCTATCCGTTGGCCCGCGATGTTTGAACCGAGGTTGATGTGCGGGCATGCTACAACTTCCAAGCTGCTTGCACTTTCGCACCATGGACGCGGTTTAATCATCACCCCAGAGACCCATGCAGAGCCGCAACAGATTGCCTTGCATGGTAACATGGGTCCGTTGCTGGCAGCTCACTGGGACGAGCGCAGCTTGATGCTTTTGGCCTCCTCCGGGGCCACGTTCCACTGCCGCGAAGCCGTTGGGTTGTGGCCTTGCCAGGACTCGCAACTCGCTCCATTGCCCTTGGGACCGGGTCCCTTCCGCGGCACCGTGGCCTTGTCACGTGCGGATGGTGACGGTGCCATGCATGCCGCGATCACCTTCCCAAAGGAACGCTCAGTGGCTATCTTTGCCCACACAGGCC
```

#### Copy Trinity 'gene' contig names

*Need whole Trinity 'gene' contig name line with length, path (rather than contig name only)*
```
enable_lmod
module load container_env perl
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" 19272_Sid_GoodCoral_Final.fasta | head
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" 19272_Sid_GoodCoral_Final.fasta > 19272_Sid_GoodCoral_Final_contigWholeLine.txt

grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" 3081_Sid_GoodSymb_Final.fasta | head
grep -P -o "(\w+_\w+\s\w+.\w+\s\w+.\[.*?[^\]]\])" 3081_Sid_GoodSymb_Final.fasta > 3081_Sid_GoodSymb_Final_contigWholeLine.txt
```

### Rename contig names to something useful
Use script fasta_name_replacer.py

*Keep record of new names and corresponding original Trinity 'gene' names*
- Create Tab delimited, 2 column table of OriginalName\tNewName (no headers).
- First column (OrignalName) is the Trinity 'gene' contig names from the above *_contig.txt file.
    - need the whole line (not only the contig Trinity name)
- Add 0's in front of digits so that the sort order is numerical. That way text sort and numerical sort are identical. 

In Excel: 

head 19272_Sid_GoodCoral_Final_contigWholeLine_renamer.txt
```
TRINITY_DN22892_c0_g1_i2 len=638 path=[0:0-445 2:446-637]	Siderastrea_Mexico_Barshis-Radice_Host_contig00001
TRINITY_DN64979_c0_g1_i6 len=578 path=[2:0-110 4:111-208 7:209-242 8:243-256 9:257-369 11:370-577]	Siderastrea_Mexico_Barshis-Radice_Host_contig00002
TRINITY_DN64831_c0_g2_i1 len=818 path=[0:0-420 2:421-817]	Siderastrea_Mexico_Barshis-Radice_Host_contig00003
TRINITY_DN484077_c0_g1_i1 len=531 path=[0:0-530]	Siderastrea_Mexico_Barshis-Radice_Host_contig00004
TRINITY_DN63583_c0_g1_i3 len=1510 path=[0:0-425 1:426-653 2:654-702 4:703-788 5:789-1509]	Siderastrea_Mexico_Barshis-Radice_Host_contig00005
TRINITY_DN91634_c0_g1_i1 len=551 path=[1:0-95 2:96-550]	Siderastrea_Mexico_Barshis-Radice_Host_contig00006
TRINITY_DN2500_c0_g1_i5 len=805 path=[0:0-76 1:77-84 4:85-250 5:251-605 7:606-626 8:627-669 10:670-751 12:752-804]	Siderastrea_Mexico_Barshis-Radice_Host_contig00007
TRINITY_DN65891_c0_g1_i1 len=755 path=[0:0-754]	Siderastrea_Mexico_Barshis-Radice_Host_contig00008
TRINITY_DN37872_c0_g1_i2 len=828 path=[0:0-583 1:584-678 3:679-827]	Siderastrea_Mexico_Barshis-Radice_Host_contig00009
TRINITY_DN36661_c0_g1_i1 len=1412 path=[1:0-1225 2:1226-1411]	Siderastrea_Mexico_Barshis-Radice_Host_contig00010
```

### Host
```
python2 fasta_name_replacer.py -i 19272_Sid_GoodCoral_Final.fasta -n 19272_Sid_GoodCoral_Final_contigWholeLine_renamer.txt -o 19272_Sid_GoodCoral_Final_renamed.fasta
```

head -4 19272_Sid_GoodCoral_Final_renamed.fasta
```
>Siderastrea_Mexico_Barshis-Radice_Host_contig00001
TTGATTTCAAGACCTGTTCACAATTGTCCCCGCAGAACCTAATCAGAAAATTTTATGGGTGAACCGGAGGTGCATAGCTTTTAAGTGTCTTGGCCAAAGAGGTTATAACATAAGGGTGAAAGTTTGAGGTGGATTCCGAACGAGGCGACCTTAAACTAAAGTTAGGCTTTTTCCTTTCTGCGTCCATCTGGAGGGCTAGACGTTCTACAAAAATTTGGTGATTTATTGGTAAAATGCTTGTCTTTTCCATTCAACCTAAATTTGAAGAATATACGATTCAGTTTACAGCCATATGTCTGAAGAACAGTCTCCCCCCGGATATCGCATTTCATGAGCAAGTGCAGGACCATCAGGCTCCTCTCCAAAGTCCACGCAGAAATCAGGATTCAACTTCAGACCACCGTTCACGGAATCTACGTCCACTTGGAGCAACATGCTGCCTTTCTTGGACATATTTGGGTAGAACTGATTGTCCCACACGCTGTACAAAGATGTGGTAACGTACAGTCGTTTTCCATCCAAACTAAGCTGAATCATTTGCGGTCCTCCTTCCACGCGCTTCCCCTTCATATAACACGGCTCAGGCTGTGCACTCAGTTCAGAGTCCTCGATCACCTTCACTGGGCCATCACGGGAT
>Siderastrea_Mexico_Barshis-Radice_Host_contig00002
AATGAGCATGTCAGTACGTAGGAGCAACGTTTCCTAGTTTTTTGTTTCTGTGCTCATTATTCTCCCTTGGGAGCAATGTAGCGCTCGTGTATTCACATTAGCTTCTTTTCCTGTCCTTTTCAGGTGTGGTTCCAAAACCGCCGTGCAAGATGGCGCAAGCGCGAAATAAAGAACAAACCTGCCCCAGTTCACCCTGCTTCGGAGAAGCGTTTGGTCAATCAAAGGGACGTAACCACCTCGAGCATGTTCCAGCCGTTGCCAAACATTTTCCCGCCACAGTACTTAGCGCCTTTGAGACCGTGGGAGACATTTTACCTGCCTTTTACGGCCAGTGCTTTGATCCCACCACAAAATATATTTTCCACTGTCATTTTAGCTTCGCATTCTTTCTTACAGTTTGTCGCCTTTAGGTGTCGTAGTGTGCTTATCGTAGCTTCCGATTGGTTCGTTTTTCCTATTTAAGATGCTGTTTTCAGTTTGTTCTCTCTCTCTCGCTTGCGACCTTAGCACTTTAGAGCTGTTTCTGCTTCGGACGTTTTAAGATTTCGTCGAGTATGCCTACTCCCGGAGACGCGGC
```

grep -c "contig" 19272_Sid_GoodCoral_Final_renamed.fasta
19272


### Symbiont
Siderastrea_Mexico_Barshis-Radice_Symbiont_contig00001 etc.

```
python2 fasta_name_replacer.py -i 3081_Sid_GoodSymb_Final.fasta -n 3081_Sid_GoodSymb_Final_contigWholeLine_renamer.txt -o 3081_Sid_GoodSymb_Final_renamed.fasta
```

head -4 3081_Sid_GoodSymb_Final_renamed.fasta
```
>Siderastrea_Mexico_Barshis-Radice_Symbiont_contig00001
GGTTGGTGTCAAGCATCATGCCCGACGTTTGCAAGCCCTAGAGGGAGTTTCAACAAGTGACATTATTTCTTCTTTGCCCCGCGAAGAGCTGGTCGCGAAGATCCGCAGCGCTTTGCCATACCTTGCCCAACTGGCAAAGGGCGGCAGTCGCCCCACTGCAGCACTTGAGAAGTCTATGATTGACACGACCACAGCGATGCCACAGGCCAAAGCCACCGCTATCCGTTGGCCCGCGATGTTTGAACCGAGGTTGATGTGCGGGCATGCTACAACTTCCAAGCTGCTTGCACTTTCGCACCATGGACGCGGTTTAATCATCACCCCAGAGACCCATGCAGAGCCGCAACAGATTGCCTTGCATGGTAACATGGGTCCGTTGCTGGCAGCTCACTGGGACGAGCGCAGCTTGATGCTTTTGGCCTCCTCCGGGGCCACGTTCCACTGCCGCGAAGCCGTTGGGTTGTGGCCTTGCCAGGACTCGCAACTCGCTCCATTGCCCTTGGGACCGGGTCCCTTCCGCGGCACCGTGGCCTTGTCACGTGCGGATGGTGACGGTGCCATGCATGCCGCGATCACCTTCCCAAAGGAACGCTCAGTGGCTATCTTTGCCCACACAGGCC
>Siderastrea_Mexico_Barshis-Radice_Symbiont_contig00002
CGCGATCGTTTGAACTTGCACCGCAAGGTTCCACGTATTGAAGCTATGGTGGTGCTCGTGGTCACCATTTTGTCCAAGCAGACCAACATTGCAGTGGCCGTGTTAGTGGGCGTCTCCATCTGCGCCATTAGCTTCGCATGGGACGCGGGCAACGACTTCAAGCTCAGCGAGTCCATGGCTGGCGACATCAAGGTTTACGATGTGGATGGTCCCTTGTTCTTCACCTCGGCCAACCGCCTCGTGAAGCTCCTGGACAGCGAGAAGGACCCACAGGAGGTGGAGCTGCGCTTCGGGGAGGCCACGCTGATGGATTTCACCTCGGTAGAGACCTTGCACAAAATCGCGCTGAACTATCAGAATGCAGGCAAAAGCATCACCTTCCACACGCTGAACCTGAGCAGCCAGAAAATCATCGAGAAGGCCAACCATTTGGTGCGCGCCATCGAGTACACCTCCTCGGATGATGCCGTGCGCTTGCCGAGCGTGCCAAGCTTCACGGAGGGCTTCCGGCACGATCCCATGACGCAGCAGGGCGTCTTCGATCGAATCGCCAGCTTCGACAGGGTCAGCAGTCCAGGCGGCAAGAAGGACAGCGCACCAAAGGATCTGGAAGTGG
```

grep -c "contig" 3081_Sid_GoodSymb_Final_renamed.fasta 
3081


*NB! Some <500 bp contigs were introduced in this process*
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 19272_Sid_GoodCoral_Final_renamed.fasta
```

```
The total number of sequences is 19272
The average sequence length is 1125
The total number of bases is 21692610
The minimum sequence length is 499
The maximum sequence length is 32630
The N50 is 1211
Median Length = 582
contigs < 150bp = 0
contigs >= 500bp = 19223
contigs >= 1000bp = 6693
contigs >= 2000bp = 1829
```

#### Filter by length threshold (> 500 bp)
```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 _500lnThresh 19272_Sid_GoodCoral_Final_renamed.fasta
```

```
Number of total seqs for 19272_Sid_GoodCoral_Final_renamed.fasta: 19272
Number of seqs over 500 for 19272_Sid_GoodCoral_Final_renamed.fasta: 19223
```

```
grep -c '>' 19272_Sid_GoodCoral_Final_renamed_500lnThresh.fasta
```
19222

**Best practice to rename with contig number**
```
cp 19272_Sid_GoodCoral_Final_renamed_500lnThresh.fasta 19222_Sid_GoodCoral_500lnThresh_Final.fasta
```

#### Check assembly
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 19222_Sid_GoodCoral_500lnThresh_Final.fasta
```

```
The total number of sequences is 19222
The average sequence length is 1127
The total number of bases is 21667620
The minimum sequence length is 500
The maximum sequence length is 32630
The N50 is 1211
Median Length = 827
contigs < 150bp = 0
contigs >= 500bp = 19222
contigs >= 1000bp = 6693
contigs >= 2000bp = 1829
```

----------------------------------------------------------------------------------------------
## Check for for overlapping contigs between de novo symbiont assembly and reference symbiont assemblies

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

De novo symbiont assembly:
3081_Sid_GoodSymb_Final_renamed.fasta

**Reciprocal Blast**
Check if any 2 contigs match same Blast result with 90% identity and whether overlapping or nearly overlapping
- 100 bp overlap
- 90% identity

clade B reference:
- 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta

clade C reference:
- 65838_davies_cladeC_feb_fixed-clean_suffixed.fasta

#### Make Blast databases

blast/2.9 (L,D)
```
enable_lmod
module load container_env blast
module avail blast
makeblastdb -in 173968_Breviolum_B5_Sid_radians_Avila-Medi_suffixed.fasta -dbtype nucl -out cladeB_Blastdb -hash_index
makeblastdb -in 65838_davies_cladeC_feb_fixed-clean_suffixed.fasta -dbtype nucl -out cladeC_Blastdb -hash_index
```

### Blast our de novo symbiont assembly as query to clade B reference

```
nano cladeB_Blast_out.sh
```

```
#!/bin/bash -l

#SBATCH -o cladeB_Blast_out.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=cladeB_Blast_out

enable_lmod
module load container_env blast
blastn -query /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/3081_Sid_GoodSymb_Final_renamed.fasta -db /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/cladeB_Blastdb -out cladeB_blastn.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
```

```
sbatch cladeB_Blast_out.sh
```

### Blast our de novo symbiont assembly as query to clade C reference

```
nano cladeC_Blast_out.sh
```

```
#!/bin/bash -l

#SBATCH -o cladeC_Blast_out.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=cladeC_Blast_out

enable_lmod
module load container_env blast
blastn -query /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref/3081_Sid_GoodSymb_Final_renamed.fasta -db /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/cladeC_Blastdb -out cladeC_blastn.outfmt5 \
        -evalue 1e-4 -num_threads 6 -max_target_seqs 1 -outfmt 5
```

```
sbatch cladeC_Blast_out.sh
```

## Parse (separate) the blast output 
Get the reads matching to the reference.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

biopython/1.75
```
enable_lmod
module load container_env python
module load biopython/1
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py cladeB_blastn_parsed.txt blastn cladeB_blastn.outfmt5
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/parse_blastnortblastx_advbioinf.py cladeC_blastn_parsed.txt blastn cladeC_blastn.outfmt5
```

grep -c "Symbiont" cladeB_blastn_parsed.txt
1721

grep -c "Symbiont" cladeC_blastn_parsed.txt
2624

### Identify "good hits" to symbiont clade B and C references
"good hits" are defined here as matching at least 90% of the read over at least 100 bp of the read

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 90per100bpmatch.txt lengthidentity 0.90 100 clade*_blastn_parsed.txt
```

```
cladeB_blastn_parsed.txt	0.90	100	Number of Good Hits:	93
cladeB_blastn_parsed.txt	0.90	100	Number of unique matches:	90
cladeC_blastn_parsed.txt	0.90	100	Number of Good Hits:	2457
cladeC_blastn_parsed.txt	0.90	100	Number of unique matches:	2152
```

**72/93 matches to clade B are also in clade C**
**Majority of our de novo symbiont assembly matches to clade C**

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/ReParseBlastbycutoffs_advbioinf_b.py 95per100bpmatch.txt lengthidentity 0.95 100 clade*_blastn_parsed.txt
```

```
cladeB_blastn_parsed.txt	0.95	100	Number of Good Hits:	10
cladeB_blastn_parsed.txt	0.95	100	Number of unique matches:	10
cladeC_blastn_parsed.txt	0.95	100	Number of Good Hits:	2407
cladeC_blastn_parsed.txt	0.95	100	Number of unique matches:	2111
```

***We decided that we will not use our symbiont de novo assembly***

***Instead will use clade B and clade C references, and analyze uniquely singly mapped reads - different DESeq analysis for each clade symbiont***


----------------------------------------------------------------------------------------------
