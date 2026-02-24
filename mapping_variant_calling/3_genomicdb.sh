#!/bin/bash

cat /home/Amel_interval.list | xargs -n 1 -P 16 -I {} bash -c \
'workspace="{}" && \
gatk GenomicsDBImport \
--java-options "-XX:ParallelGCThreads=2" \
--genomicsdb-workspace-path $workspace \
--sample-name-map sample.sample_map \
--batch-size 90 \
--intervals {} && \'
