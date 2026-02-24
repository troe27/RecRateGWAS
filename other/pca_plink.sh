#!/bin/bash -l

in_file="concat_chr_GQ20_Fmissing_hetALL.vcf.gz"
out_file=$in_file"_plink"

plink --vcf $in_file --allow-extra-chr --const-fid --vcf-half-call m --recode --out $out_file

in_file2=$out_file
out_file2=$in_file2"_pca"

plink --file $in_file2 --allow-extra-chr --pca tabs --out $out_file2

