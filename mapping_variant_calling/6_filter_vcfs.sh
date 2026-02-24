#!/bin/bash

# Define the process you want to monitor
PROCESS_NAME="setGT"  # Replace with the process name or ID
PROCESS_PID=$(pgrep -f "$PROCESS_NAME")  # Get the PID of the process

# Function to check if the process is still running
is_process_running() {
    if ps -p $PROCESS_PID > /dev/null; then
        return 0  # Process is still running
    else
        return 1  # Process is not running
    fi
}

# Check if the process is running, wait for it to finish
echo "Monitoring the process: $PROCESS_NAME (PID: $PROCESS_PID)"
while is_process_running; do
    echo "Waiting for $PROCESS_NAME to be done"
    sleep 300  # Check every 5 minutes
done

# Process finished, now run the commands sequentially
echo "Process finished. Running other commands."

echo "All commands have completed."

# Run GATK VariantFiltration
    gatk VariantFiltration -R /home/drone_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
        -V ${file_id}_drone_only_biallelicSNPs.vcf.gz \
        --filter-name "QD" --filter-expression "QD < 20.0" \
        --filter-name "FS" --filter-expression "FS > 60.0" \
        --filter-name "MQ" --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum" --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
        --filter-name "SOR" --filter-expression "SOR > 3.0" \
        -O ${file_id}.SNP.filtered.vcf.gz

    tabix -f ${file_id}.SNP.filtered.vcf.gz

	# recalculate GQ-values in order to reflect that the sequences are actually haploid
    python3 /home/gvcf/1185results/hapGQ.py ${file_id}.SNP.filtered.vcf.gz ${file_id}_GQ.vcf.gz
    tabix -f ${file_id}_GQ.vcf.gz

    bcftools +setGT ${file_id}_GQ.vcf.gz -- -t q -n . -i "GQ<20" | bgzip > ${file_id}_GQ20.vcf.gz
    tabix -f ${file_id}_GQ20.vcf.gz
    
    bcftools filter -i "INFO/DP < 20000" /home/gvcf/1185results/${file_id}_GQ20.vcf.gz -O z9 -o /home/gvcf/1185results/${file_id}_GQ20_DP.vcf.gz; vcftools --gzvcf /home/gvcf/1185results/${file_id}_GQ20_DP.vcf.gz --max-missing 0.9 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout | bcftools view -i "COUNT(GT==\"het\")<=10" -Ou | bcftools +setGT -- -t q -n . -i "GT==\"het\"" | bgzip > /home/gvcf/1185results/yapp_1224samples/${file_id}_GQ20_Fmissing_hetALL.vcf.gz

