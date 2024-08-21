#!/bin/bash
#SBATCH --account=ELLIS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --job-name=bulk_cu_nxnxn
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=12:00:00

. /etc/profile.d/modules.sh


module purge
module load rhel8/default-icl castep

mpirun castep bulk_cu_nxnxn