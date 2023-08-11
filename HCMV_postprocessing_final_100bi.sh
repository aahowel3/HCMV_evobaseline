#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 2-00:00:00   # time in d-hh:mm:ss
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

for i in {1..25}
do
for j in {1..4}
do
cd replicate_${i}_0${j}

for file in congenital_plasma*.100.output.ms
do
base=$(basename $file .100.output.ms)
strip="${base#*.}"

#echo replicate_${i}_${j}
Rscript ../HCMV_postprocessing_final_100bi.R $file congenital_plasma.${strip}.1000.output.ms congenital_urine.${strip}.100.output.ms congenital_urine.${strip}.1000.output.ms  >> 100bicombos_2.csv
#Rscript HCMV_postprocessing_final.R 
done
cd ../
done
done
