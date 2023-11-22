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

for i in {1..10}
do
#mkdir replicate_"$i"
cd replicate_"$i"
        #generate all .ms and .fix files according to the list of param combinati9ons
        while IFS="," read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11
                do
#	echo $col1 col2 col3 col4 col5 $col6 $col7 $col8 $col9 $col10 $col11
        coll1=$(sed -e 's/^"//' -e 's/"$//' <<<"$col1")
#	coll11=$(sed -e 's/^"//' -e 's/"$//' <<<"$col11")
	coll7=$(sed 's/./&0/5' <<<"$col7") 
	repls=2.0
	coll9="${col9/2/$repls}"
	colll9=$(sed 's/./&0/5' <<<"$coll9")
	coll10=$(sed -e 's/[0]*$//g' <<<"$col10")
	coll11=$(sed -e 's/[0]$//g' <<<"$col11")
	fname="congenital."${coll1}".${colll9}.${coll7}.${coll10}.${coll11}.fix.txt"
#	echo $fname
	if [[ -f "$fname" ]]
        then
                continue
        fi 
#        echo "missing"
		slim -d DFE="$col1" -d f0="$col2" -d f1="$col3" -d f2="$col4" -d f3="$col5" -d recombrate="$col6" -d recomb="$col7" -d murate="$col8" -d mu="$col9" -d progeny="$col10" -d gr="$col11" ../HCMV_burnin_final.slim
        done < ../params_burnin_nor.txt
	cd ../
done
