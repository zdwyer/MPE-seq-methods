---
title: "High-precision detection of pre-mRNA splicing isoforms using Multiplexed Primer Extension sequencing"
author: "Zachary W Dwyer, Michael A Gildea, and Jeffrey A Pleiss"
---

This document contains the data processing and analysis done for "High-precision detection of pre-mRNA splicing isoforms using Multiplexed Primer Extension sequencing" published in Methods. All custom scripts can be found in the "scripts" folder of this project. Questions and comments can be sent to Zach Dwyer at zwd2@cornell.edu.

# Figure 1

## Data Processing

### *Schizosaccharomyces pombe*
The 2.30 release of the *Schizosaccharomyces pombe* genome was downloaded from [Pombase](ftp://ftp.pombase.org/pombe/genome_sequence_and_features/).

## Data Analysis

# Figure 2

## Data Processing

### Genome and Annotation Files
The R64-2-1 release of the *Saccharomyces cerevisiae* genome was downloaded from the [Saccharomyces Genome Database](https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/). The genome sequence was manually removed from the feature file and saved as sc_feature_R64-2-1.gff. The chromosome names of the genome file were manually renamed to match the feature file and saved as sc_genome_R64-2-1.fa.

#### Build HISAT index:
```
python sc_extract_exons_for_hisat.py sc_feature_R64-2-1.gff > sc_exons.txt
python sc_extract_introns_for_hisat.py sc_feature_R64-2-1.gff > sc_introns.txt

hisat2-build --ss sc_introns.txt --exon sc_exons.txt sc_genome_R64-2-1.fa sc_index_R64-2-1
```

### Sample Summary

|**Source** | **Organism** | **Strain** | **Replicate** | **Library Preparation** | **Read Length** | **Sample Name** | **Reads**   |
|-----------|--------------|------------|:-------------:|:-----------------------:|-----------------|-----------------|------------:|
|[SRR7208763](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7208763) | *S. cerevisiae* | BY4741     | A             | RNA-seq                 | 101             | sc_RNA_A        | 27,248,124  | 
|[SRR7208770](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7208770) | *S. cerevisiae* | BY4741     | B             | RNA-seq                 | 101             | sc_RNA_B        | 30,116,992  |
|[ID]()         | *S. cerevisiae* | BY4741     | A             | MPE-seq                 | 60+15           | sc_MPE_A        | 4,986,024   | 
|[ID]()         | *S. cerevisiae* | BY4741     | B             | MPE-seq                 | 60+15           | sc_MPE_B        | 5,803,915   |

### Downsample to 5 Million Reads
```

```

### Remove PCR Duplicates [MPE-seq Only]
```
python compress_UMI.py -n 7 -1 sc_MPE_A_R1.fastq.gz -2 sc_MPE_A_R2.fastq.gz > sc_MPE_A.fastq.gz
```

### Trimming
Reads are trimmed with (Trimmomatic version 0.35) to remove any read through into sequencing adapters which will interfear with alignment. RNA-seq libraries were prepared with TruSeq adapters while MPE-seq libraries were prepared with Nextera adapters. MPE-seq libraries were required to be at least 26 bases long post trimming to remove any reads that come from unextended primers during reverse transcription. TruSeq3-SE.fa and NexteraPE-PE.fa came with Trimmomatic and can additionally be found in the "resources" folder of the project. Example commands:

#### RNA-seq
```
java -jar trimmomatic.jar SE -phred33 sc_RNA_ds_A.fastq.gz sc_RNA_A_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:10
```
#### MPE-seq
```
java -jar trimmomatic.jar SE -phred33 sc_MPE_ds_A.fastq.gz sc_MPE_A_trimmed.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 MINLEN:26
```
Summary of trimming:

| **Sample Name** | **Strain** | **Replicate** | **Library Preparation** | **Reads**  | **Survived** | **Dropped** | **Percent Survived** |
|-----------------|------------|:-------------:|-------------------------|-----------:|-------------:|------------:|---------------------:|
| sc_RNA_A        | BY4741     | A             | RNA-seq                 | 5,000,000  | 4,988,092    | 11,908      | 99.76                |
| sc_RNA_B        | BY4741     | B             | RNA-seq                 | 5,000,000  | 4,988,824    | 11,176      | 99.78                |
| sc_MPE_A        | BY4741     | A             | MPE_seq                 | 3,645,443  | 3,642,108    | 3,335       | 99.91                |
| sc_MPE_B        | BY4741     | B             | MPE_seq                 | 3,467,732  | 3,465,210    | 2,522       | 99.93                |

### Alignment
Reads were aligned with HISAT version 2.1.0. Reads with a mapping quality of less than 5 were filtered out (they typically arise from multimapping reads). Example comand:
```
hisat2 --max-intronlen 2000 --summary-file sc_RNA_A_AlignmentSummary.txt --new-summary --no-unal -p 4 -x sc_index_R64-2-1 -U sc_RNA_A_trimmed.fastq.gz | samtools view -bh -q 5 - | samtools sort - -o sc_RNA_A.bam
```
Summary of alignment:

|**Sample Name** | **Strain** | **Replicate** | **Library Preparation** | **Reads**  | **Unaligned** | **1 Alignment** | **Multiple Alignments** |
|----------------|------------|:-------------:|-------------------------|-----------:|--------------:|----------------:|------------------------:|
| sc_RNA_A       | BY4741     | A             | RNA-seq                 | 4,988,092  | 516,491       | 4,205,157       | 266,444                 |
| sc_RNA_B       | BY4741     | B             | RNA-seq                 | 4,988,824  | 396,795       | 4,309,588       | 282,441                 |
| sc_MPE_A       | BY4741     | A             | MPE-seq                 |   |       |       |                 |
| sc_MPE_B       | BY4741     | B             | MPE-seq                 |   |       |       |                 |

### Feature Counting
Mature alignments and premature reads were counted with a custom script that is based off of HTSeq-count. Mature alignments are those that cross an exon-exon junction. Premature alignments are either across an exon-intron or intron-exon boundary or are completely intronic. Results can be found in the "data" folder of this project.
```
python feature_count2.py -i sc_RNA_A.bam -f sc_intron_ranges.bed > sc_RNA_A_counts.txt
```

## Analysis

