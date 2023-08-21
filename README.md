 patientdata.sh aligns and does variant calling on B103 datasets
 HCMV_filteringB103.R gets the number of segregating sites of B103 6mo plasma and kidney compartment samples at 100X and 1000x subsampling 
 HCMV_postprocessing R and sh pairs then evaluate the number of segregating sites in the simulations against the empirical number of segergating sites in B103 on filter thresholds of 100 or 10 for bialleic and triallelic frequencies

 get the sumstatcombos.txt file to iterate through only viable combinations per replicate to perform SC2 on - 
 for i in replicate*; do cat $i/100bi_10tricombos_2.csv | grep "plasma_100_bi" | grep "good" | awk '{print $2}' > $i/sumstatcombos.txt; done 
 
 HCMV_convert_to_ms_notris.R and HCMV_summarystats.sh then create new 2% MAF filtered .ms files from the original simulation .ms files and run Terbot's sc2 summary statistics scripts on the filtered .ms files


In directory: /scratch/aahowel3/simulations_replicates/good_combos_sumstats] 
concatfor_summarystats_plotting.sh - takes all the .filter.csv from each replicate and combines them per combination
summarystats_plotting.sh - loops each param combination through plot_summarystats_updated_updated.R to get a 5 page 4 plot per page 1 document per parameter combination file of all 4 compartments across the 5 summary statistics 
