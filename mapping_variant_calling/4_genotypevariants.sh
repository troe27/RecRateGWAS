#!/bin/bash

PROCESS_NAME="GenomicsDBImport"  # Replace with the process name or ID
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

cat /home/gvcf/regions.list | xargs -n 1 -P 16 -I {} bash -c '

    chrom_name=$(echo {} | cut -d ":" -f 1)
    
    gatk GenotypeGVCFs \
        -R /home/drone_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
        -V gendb://$chrom_name \
        -L {} \
        --allow-old-rms-mapping-quality-annotation-data \
        -O /home/gvcf/1185results/chromosome_fragments_g_vcf/{}.vcf.gz
'
