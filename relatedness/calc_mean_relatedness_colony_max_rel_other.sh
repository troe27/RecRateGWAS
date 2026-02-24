#!/bin/bash -l

infile="concat_chr_GQ_relatedness_set4.relatedness"
outfile=$infile"_mean"

>$outfile
printf "colony\tindv\tmean_rel_col\tmax_rel_other\tmax_rel_other_indv1\tmax_rel_other_indv2\n" >$outfile

indv_set1=`cat $infile | tail -n +2 | awk '{print $1}' | sort | uniq`

for indv in $indv_set1; do
	colony=`echo $indv | awk -F "_" '{print $1"_"$2}'`
	mean_rel_col=`cat $infile | grep -w $indv | grep $colony"_.*"$colony"_" | awk 'BEGIN {rel=0; i=0} {rel+=$3;i++} END {print rel/i}'`
	printf "%s\t%s\t%f\t" $colony $indv $mean_rel_col >> $outfile
	cat $infile | grep -w $indv | grep -v $colony"_.*"$colony"_" | awk 'BEGIN {max_rel=-20;max_indv1="NA";max_indv2="NA"} {if ($3>max_rel) {max_rel=$3;max_indv1=$1;max_indv2=$2} } END {print max_rel"\t"max_indv1"\t"max_indv2}' >> $outfile
done


