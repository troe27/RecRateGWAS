#!/bin/bash -l

# Run vcftools --relatedness in order to produce input files

infile="concat_chr_GQ_relatedness_set4.relatedness"
outfile=$infile"_matrix"
>$outfile

printf "indv_vs_indv\t" > $outfile

indv_0=`cat $infile | tail -n +2 | head -n 1 | awk '{print $1}'`

cat $infile | grep -w "^"$indv_0 | awk '{print $2}' | tr '\n' '\t' >> $outfile
printf "\n" >> $outfile

dummy=""

for indv in `cat $infile | tail -n +2 | awk '{print $1}' | uniq`; do
	rel_data=`cat $infile | grep -w "^"$indv | awk '{print $3}' | tr '\n' '\t'`
	printf "%s\t%s%s\n" "$indv" "$dummy" "$rel_data" >> $outfile
	dummy=`printf "%s0\t" "$dummy"`
done


