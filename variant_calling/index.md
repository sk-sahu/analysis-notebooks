# Variant Calling
Starting from raw reads (.fastq) to calling variants (.vcf) and their downstream analysis.

## Tools Required:
* BWA 
* Picard 
* GATK-4 
* SnpEFF 
* Samtools 

## Tool installation
Install the required tools in a conda environment using “variant_call.yml” file.

But before that, you need to have conda installed in the system. 
If it's already not installed. To install conda, from the command line:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda3
```

Once the conda is installed. Now let's install the tools. This will take time because all the tools will be installed at once.

```
conda env create -f environment.yml
```

Now, If you type `conda env list` you should be able to see an environment called - variant_call. Change to that environment by entering

```
!conda activate variant_call
```

## Start Analysis

### Create Reference Genome Index

```
bwa index ref_genome.fa
samtools faidx ref_genome.fa
picard CreateSequenceDictionary R= ref_genome.fa O= ref_genome.fa.dict 
```

### Step1: Alignment – Map to Reference

```
bwa mem -M -t 16 ref_genome.fa read1.fq read2.fq > aln.sam
```

### Step2: Generate sorted BAM file

```
picard SortSam INPUT=aln.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate
```

### Step3: Check Alignment Summary

```
picard CollectAlignmentSummaryMetrics R=ref_gneome.fa I=sorted_reads.bam O=alignment_metrics.txt

picard CollectInsertSizeMetrics INPUT=sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf

samtools depth -a sorted_reads.bam > depth_out.txt
```

### Step4: Mark Duplicates 

```
picard MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt
```

### Step:5 Add Read Groups

```
picard AddOrReplaceReadGroups I=dedup_reads.bam O=final.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
```

### Step6: Build BAM Index

```
picard BuildBamIndex INPUT=final.bam
```

### Step7: Call Variants

```
gatk HaplotypeCaller -R ref_genome.fa -I final.bam -O raw_variants.vcf
```

### Step8: Extract SNPs & Indels

#### For SNP

```
gatk SelectVariants -R ref_genome.fa -V raw_variants.vcf -select-type SNP -O raw_snps.vcf
```

#### For Indel

```
gatk SelectVariants -R ref_genome.fa -V raw_variants.vcf -select-type INDEL -O raw_indels.vcf
```
