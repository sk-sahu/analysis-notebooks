# Here we will demonstarte variant calling using freebayes
The steps-1 to 6 will be the same as of mentioned [here](https://github.com/sk-sahu/notebooks/blob/master/vartiant_call.md)

## Enviroment
This suppose your are in the `variant_call` conda enviroment. Otheriwise activate it - `conda activate variant_call`.
Then install freebayes in this enviroment.
```
conda install freebayes==1.3.1 -c bioconda
```

## Step7: Call Variants
```
freebayes -f ref_genome.fa final.bam > raw_variants.vcf
```

## Filter VCFs
```
vcftools --vcf raw_variants.vcf --minQ 20 --recode --recode-INFO-all --out filtred_variants.vcf
```

## Annotation with known variants
Download required refernece dbsnp database in vcf format (for human from  here: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/)

Compress and index your `raw_variants.vcf` file to be used by next step.
```
bgzip filtred_variants.vcf 
tabix filtred_variants.vcf.gz
```
Also index downloaded `reference.vcf`
```
tabix reference.vcf.gz
```
Annotate using `bcftools`
```
bcftools annotate -c ID \
        -a reference.vcf.gz filtred_variants.vcf.gz \
        > filtred_variants_annot.vcf
```


# Reference
[Variant-Calling-using-freebayes-and-Annotation](https://github.com/CBC-UCONN/Variant-Calling-using-freebayes-and-Annotation)