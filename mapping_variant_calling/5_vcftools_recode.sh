#!/bin/bash

cat /home/Amel_interval.list | xargs -n 1 -P 16 -I {} bash -c 'vcftools --gzvcf {}.SNP.filtered.vcf.gz --remove-indels --remove-filtered-all --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --stdout | gzip -c > {}.remove.filtered.vcf.gz'
