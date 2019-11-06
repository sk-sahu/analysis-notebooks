# Compare VCF files

## Install the required tools
It supposed that you already have `variant_call` conda enviroment in your system from the previous analysis.
```
conda activate variant_call
conda install -c tabix bioconda bcftools bedtools vcftools perl-vcftools-vcf
```

## Explore vcf file
First of all comopress all the vcf files and index them (ignore if they are already done)
```
bgzip gatk_snps.vcf
tabix gatk_snps.vcf.gz
```
### How many snps vcf file have
```
bcftools view -v snps gatk_snps.vt.vcf | grep -v "^#" | wc -l
```

### Look for number of common snps in both the file
```
bedtools intersect -u -a gatk_snps.vcf.gz -b freebayes_snps.vcf.gz | wc -l
```

### calculate Jaccard index
Jaccard score indicates the operlap of snp cordinates.
```
bedtools jaccard -a gatk_snp.vcf.gz -b freebayes_snp.vcf.gz
```

### Using vcf-compare
this gives nice summary of all the common and uniue snps belongs to each vcf file.
```
vcf-compare freebayes_snp.vcf.gz gatk_snp.vcf.gz
```
Take the given numbers in `VN` and make a venn diagram.

### Using bcftools
```
bcftools isec bcftools_snps.vcf.gz gatk_snps.vcf.gz -p isec_dir
```

