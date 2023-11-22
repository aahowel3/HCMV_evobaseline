#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 4-00:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user=%u@asu.edu # Mail-to address
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module load mamba
source activate HCMV
module load r-4.2.2-gcc-11.2.0

for file in ../F2_rescaled/congenital_plasma*DFE1*1.0.with*.100.output.ms
do
base=$(basename $file .100.output.ms)
strip="${base#*.}"
Rscript HCMV_segsites_sims2.R $file ../F2_rescaled/congenital_plasma.${strip}.1000.output.ms \
../F2_rescaled/congenital_urine.${strip}.100.output.ms ../F2_rescaled/congenital_urine.${strip}.1000.output.ms  >> combos1.csv
done

