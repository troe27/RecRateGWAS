#!/bin/bash -l

# Main directories as variables
# SEQDIR assumed to contain all fastq-files

# Output: 	Sorted, indexed BAM-file for each sample
# 		BAM-stats for all samples written to the same file (bamtools_stats.txt in SEQDIR)

# Update this with the full path to the directory where you have the sample fastq files
SEQDIR="/home/drone_data/fastq"
OUTDIR="/home/drone_data/bam"

# Reference sequence
# Update this with the full path to the directory where you have the ref genome
REFDIR="/home/drone_data/ref"
# Update this with the name of the ref fasta file
REF=$REFDIR"/GCF_003254395.2_Amel_HAv3.1_genomic.fna"

map_reads () {

	# Load the bwa module to make bwa available
	# Load samtools, bamtools

	module load bioinfo-tools
	module load bwa
	module load samtools
	module load bamtools
	
	# Enter directory with data...
	
	cd $SEQDIR/
	
	# Create file for bamtools stats
	STAT_FILE=$OUTDIR/"bamtools_stats.txt"
	>$STAT_FILE
	
	for SAMP_DIR in `cat samp_list.txt`; do 

		    			echo "MAPPING $SAMP_DIR ..."

                    			# Run bwa across all allocated cores with the new mem
                    			# algorithm, and pipe the output to samtools to make a
                    			# sorted BAM file
                    
                    			time bwa mem -t 12 $REF \
						${SAMP_DIR}*"_R1"*"_001.fastq.gz" \
						${SAMP_DIR}*"_R2"*"_001.fastq.gz" \
						2> ${OUTDIR}/${SAMP_DIR}.bam.bwa.log \
						| samtools view -b | samtools sort -o ${OUTDIR}/${SAMP_DIR}.sorted.bam -O BAM
                    
                    			# Index the BAM file
                   			cd $OUTDIR 
                    			samtools index ${SAMP_DIR}.sorted.bam

				# If sorted BAM-file was created in previous step, or existed from before: extract the bamstats and append to file
				if [ -e ${SAMP_DIR}.sorted.bam ]; then

					echo "Running bamtools stats for ${SAMP_DIR}.sorted.bam"
					echo "${SAMP_DIR}" >> $STAT_FILE
					bamtools stats -in ${SAMP_DIR}.sorted.bam >> $STAT_FILE

				else

					echo "BAM-file ${SAMP_DIR}.sorted.bam not found. No stats produced."

				fi
        
				cd $SEQDIR
	done

}

# Run the function map_reads above 

map_reads






