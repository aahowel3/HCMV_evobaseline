 patientdata.sh aligns and does variant calling on B103 datasets
 HCMV_filteringB103.R gets the number of segregating sites of B103 6mo plasma and kidney compartment samples at 100X and 1000x subsampling 
 HCMV_postprocessing R and sh pairs then evaluate the number of segregating sites in the simulations against the empirical number of segergating sites in B103 on filter thresholds of 100 or 10 for bialleic and triallelic frequencies
 HCMV_convert_to_ms.R and HCMV_summarystats.sh then create new 2% MAF filtered .ms files from the original simulation .ms files and run Terbot's sc2 summary statistics scripts on the filtered .ms files


 get the sumstatcombos.txt file to iterate through only viable combinations per replicate to perform SC2 on - loop through per replicate
 cat 100bi_10tricombos.csv | grep "plasma_100_bi" | grep "good" | awk '{print $2}'
 
