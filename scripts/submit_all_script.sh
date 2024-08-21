#!/bin/bash

# Loop from 1 to 10
for i in {1..10}; do
    # Navigate to the directory
    cd "./bulk_cu_${i}x${i}x${i}" || { echo "Failed to cd into directory"; exit 1; }

    # Submit the SLURM job
    sbatch submit_job.sh

    # Return to the previous directory
    cd .. || { echo "Failed to cd back"; exit 1; }

    # Print the status message
    echo "${i} done"
done

