#!/bin/bash

# Define the process you want to monitor
PROCESS_NAME="bash -c"  # Replace with the process name or ID
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

cat /home/hieu/try_hieu/Amel_interval.list | xargs -n 1 -P 16 -I {} bash -c '
    file_id={}
    tabix /home/drone_data/gvcf/1185results/yapp_1224samples/${file_id}_GQ20_Fmissing_hetALL.vcf.gz
    yapp sperm --err 0.01 --rho 26 ${file_id}_GQ20_Fmissing_hetALL

echo "All commands have completed."
'
