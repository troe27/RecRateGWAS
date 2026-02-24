#!/bin/bash -l

# Mark duplicate reads and add read groups to samples
# Also rename bam-files to names that better represents the samples (eg bee colony id)
# ***************Switched order of addreplacereadgroups and markduplicates*************

# Final output: new_name.sorted.marked_dupes.rg.bam, ...bam.bai

# The addition of read groups probably needs to be updated depending on the actual sequencing data structure
# In this example: one sample per library, all on the same plate, same lane, same flowcell

#Path to the directory where you have the bam-files
SEQDIR="/home/drone_data/bam"

#METADATA is a file that contains the following columns for each sample, separated by a tab character
#lane; date; flowcell; old sample name; library name; new name (per individual)
#there is one row for each sample/library

METADATA=$SEQDIR/metadata_new.txt

num_samples=`cat $METADATA | wc -l`

#Store information from the METADATA file as arrays that can be acessed in the for-loop below

lanes_all=`awk '{  print $1 }' $METADATA`
dates_all=`awk '{  print $2 }' $METADATA`
flowcells_all=`awk '{  print $3 }' $METADATA`
samp_old_all=`awk '{  print $4 }' $METADATA`
library_all=`awk '{  print $5 }' $METADATA`
indv_all=`awk '{  print $6 }' $METADATA`

declare -a lanes_array
declare -a dates_array
declare -a flowcells_array
declare -a samp_old_array
declare -a library_array
declare -a indv_array

j=0
for lane in ${lanes_all[*]}; do
	lanes_array[j]=$lane
	j=$(( $j + 1 ))
done

j=0
for seq_date in ${dates_all[*]}; do
	dates_array[j]=$seq_date
	j=$(( $j + 1 ))
done

j=0
for flowcell in ${flowcells_all[*]}; do
	flowcells_array[j]=$flowcell
	j=$(( $j + 1 ))
done

j=0
for samp_old in ${samp_old_all[*]}; do
	samp_old_array[j]=$samp_old
	j=$(( $j + 1 ))
done

j=0
for library in ${library_all[*]}; do
	library_array[j]=$library
	j=$(( $j + 1 ))
done

j=0
for indv in ${indv_all[*]}; do
	indv_array[j]=$indv
	j=$(( $j + 1 ))
done


function arrange_read_groups {

	module load bioinfo-tools
	module load picard/1.118
	module load samtools

	# Parameters for specifying read groups 

	# RGID = Identifier (same run and library)  unique for each run (e.g. has lane)
	# RGPL = Platform (Solid, Illumina etc)     same across project
	# RGPU = Platform unit (run barcode)        ~ unique for each library (barcode)
	# RGLB = Library                            unique for each library (shared among lanes)
	# RGSM = Sample                             unique for each library (shared among lanes)

	SORT_ORDER=coordinate

	cd $SEQDIR
	
	for sample in `seq 0 $(( $num_samples - 1 ))`; do

		old_name=${samp_old_array[sample]}

		library=${library_array[sample]}

		indv_name=${indv_array[sample]}

		lane=${lanes_array[sample]}

		seq_date=${dates_array[sample]}

		flowcell=${flowcells_array[sample]}

		new_name=$indv_name
            
            	RGFIXED=""
                
                # --------------
                # Add readgroups
		# --------------
                
		echo "Adding readgroups for $old_name.sorted.bam ..."

                RGID=$seq_date""$flowcell""$lane
        	RGPL="ILLUMINA"
                RGPU=$flowcell"."$lane"."$library
                RGLB=$library
                RGSM=$indv_name
                
                time java -Xmx8g -jar "/proj/programs/picard.jar" \
                	AddOrReplaceReadGroups \
                    	INPUT=$old_name.sorted.bam \
                    	OUTPUT=$new_name.sorted.rg.bam \
                    	SORT_ORDER=$SORT_ORDER \
                        RGID=$RGID \
                        RGPL=$RGPL \
                        RGPU=$RGPU \
                        RGLB=$RGLB \
                        RGSM=$RGSM

                    
                
                # ---------------
                # Mark duplicates
                # ---------------
                
		INPUT=$new_name".sorted.rg.bam"

                echo "Marking duplicates for $INPUT ..."
                
                time java -Xmx8g -jar "/proj/programs/picard.jar" MarkDuplicates \
			INPUT=$INPUT OUTPUT=$new_name.sorted.marked_dupes.rg.bam METRICS_FILE=$new_name.sorted.bam.marked_dupes.rg.metrics.csv

                # Index the file

                samtools index $new_name.sorted.marked_dupes.rg.bam
                samtools idxstats $new_name.sorted.marked_dupes.rg.bam > $new_name.sorted.marked_dupes.rg.bam.idxstats.csv
    
    done	

}

arrange_read_groups
	
