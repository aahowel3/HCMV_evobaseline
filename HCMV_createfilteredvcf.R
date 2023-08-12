####################
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/patient_data/good_variants_final/")

library('VariantAnnotation')
#convert vcf dataframe from B103 into filtered output vcf to use in sc2 summarystats script
vcf <- readVcf("../475_sam_dupl_rm_no5_plasma6mo.vcf")
vcf <- readVcf("../475_sam_dupl_rm_no5_urine6mo.vcf")

df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3 <- merge(df1,df2, by="row.names")


#SPLIT AO
for(i in 1:length(df3$AO)){
  df3$AO[i] <- paste(df3$AO[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$AO=gsub( "\\c", "", as.character(df3$AO))
df3$AO=gsub( "\\(", "", as.character(df3$AO))
df3$AO=gsub( "\\)", "", as.character(df3$AO))
df3$AO=gsub( " ", "", as.character(df3$AO))
df3$AO=gsub( ":", ",", as.character(df3$AO))

#need missing entries to be NA stringploit was just putting first value in third col??? 
df3 = df3 %>% 
  separate(col = AO, sep=",",into= str_c('AO',1:3 )) # color1 through color7 column names


#SPLIT SAR
for(i in 1:length(df3$SAR)){
  df3$SAR[i] <- paste(df3$SAR[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$SAR=gsub( "\\c", "", as.character(df3$SAR))
df3$SAR=gsub( "\\(", "", as.character(df3$SAR))
df3$SAR=gsub( "\\)", "", as.character(df3$SAR))
df3$SAR=gsub( " ", "", as.character(df3$SAR))
df3$SAR=gsub( ":", ",", as.character(df3$SAR))

df3 = df3 %>% 
  separate(col = SAR, sep=",",into= str_c('SAR',1:3 )) # color1 through color7 column names

#SPLIT SAF
for(i in 1:length(df3$SAF)){
  df3$SAF[i] <- paste(df3$SAF[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$SAF=gsub( "\\c", "", as.character(df3$SAF))
df3$SAF=gsub( "\\(", "", as.character(df3$SAF))
df3$SAF=gsub( "\\)", "", as.character(df3$SAF))
df3$SAF=gsub( " ", "", as.character(df3$SAF))
df3$SAF=gsub( ":", ",", as.character(df3$SAF))

df3 = df3 %>% 
  separate(col = SAF, sep=",",into= str_c('SAF',1:3 )) # color1 through color7 column names


#SPLIT RPR reads placed right
for(i in 1:length(df3$RPR)){
  df3$RPR[i] <- paste(df3$RPR[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$RPR=gsub( "\\c", "", as.character(df3$RPR))
df3$RPR=gsub( "\\(", "", as.character(df3$RPR))
df3$RPR=gsub( "\\)", "", as.character(df3$RPR))
df3$RPR=gsub( " ", "", as.character(df3$RPR))
df3$RPR=gsub( ":", ",", as.character(df3$RPR))

df3 = df3 %>% 
  separate(col = RPR, sep=",",into= str_c('RPR',1:3 )) # color1 through color7 column names

#SPLIT RPRL reads placed left 
for(i in 1:length(df3$RPL)){
  df3$RPL[i] <- paste(df3$RPL[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$RPL=gsub( "\\c", "", as.character(df3$RPL))
df3$RPL=gsub( "\\(", "", as.character(df3$RPL))
df3$RPL=gsub( "\\)", "", as.character(df3$RPL))
df3$RPL=gsub( " ", "", as.character(df3$RPL))
df3$RPL=gsub( ":", ",", as.character(df3$RPL))

df3 = df3 %>% 
  separate(col = RPL, sep=",",into= str_c('RPL',1:3 )) # color1 through color7 column names



#so using df3 as the true number of visible sites - how many of them are polymorphic? 
#this is where we can then use qual etc 
df=df3[df3$DP > 100,]
df=df[df$QUAL > 0,]



#cut down to 100 frequency
df$RO_down <- round(100*(df$RO/df$DP))
df$AO1_down <- round(100*(as.numeric(df$AO1)/df$DP))
df$AO2_down <- round(100*(as.numeric(df$AO2)/df$DP))
df$AO3_down <- round(100*(as.numeric(df$AO3)/df$DP))

#get rid of sites where ref is greater than 98 - if 99 other sites cant be greater than 2%
df=df[df$RO_down <= 98,]
#get rid of sites where no alt is at least greater than 2% but also less than 99
df=df[df$AO1_down <=98 & df$AO1_down > 1 | df$AO2_down <=98 & df$AO2_down > 1 | df$AO3_down <=98 & df$AO3_down > 1,] 
df=df[complete.cases(df$Row.names),]

##THIS - THIS IS THE RIGHT IDEA 
#NEEDS TO BE AT LEAST TWO SITES MEETTHING THESE CONDITIONS - NOT OROROR 
#OTHERWISE ITS NOT TRULY POLYMORPHIC 
#SOMETHINGS BEING WEIRD WITH THE NA THO
df=as.data.frame(df, stringsAsFactors = FALSE)
#perfect! just throws a warning not an error!!
df[is.na(df)] <- 0

###whoops whoops cant od this 
#all those sites still need to be supported by at least one read one ither sites and at that site 
##tris=subset(df5,  rowSums(df5[65:68] >= 2) == 3)
#PREV CHUNK WAS FOR TOTAL THIS NOW BREAKS DOWN BY BI/TRI/QUAD ALLELEIC
# == INSTEAD OF >=
df5 <- df[rowSums(cbind((df$RO_down >= 2 & df$SRR > 0 & df$SRF > 0), (df$AO1_down >= 2 & df$SAR1 > 0 & df$SAF1 > 0), 
                        (df$AO2_down >= 2 & df$SAR2 > 0 & df$SAF2 > 0),  
                        (df$AO3_down >= 2 & df$SAR3 > 0 & df$SAF3 > 0))) == 2,] #CHANGE HERE FOR BIALLEICS (2) #CHANGE HERE FOR TRIALLEICS (3) #CHANGE HERE FOR QUADALLEICS (4) 
df5=df5[complete.cases(df5$Row.names),]
df5=subset(df5, ((RO_down >= 2) & rowSums(cbind((df5$AO1_down >= 2 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                                (df5$AO2_down >= 2 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                                (df5$AO3_down >= 2 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 1) | #CHANGE HERE FOR BIALLEICS (1) #CHANGE HERE FOR TRIALLEICS (2) #CHANGE HERE FOR QUADALLEICS (3) 
             ((RO_down < 2) & rowSums(cbind((df5$AO1_down >= 2 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                            (df5$AO2_down >= 2 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                            (df5$AO3_down >= 2 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 2)) #CHANGE HERE FOR BIALLEICS (2) #CHANGE HERE FOR TRIALLEICS (3) #CHANGE HERE FOR QUADALLEICS (4) 
df5=df5[complete.cases(df5$Row.names),]

#pick the largest of the AOs as the alt
#100-that number becomes RO
#should account for if there's multiple AOs as highest
library(dplyr)
df5=df5 %>% 
  mutate(AO_down = pmax(AO1_down, AO2_down, AO3_down))

df5$RO_downf=100-df5$AO_down

#in new vcf - replace RO values here with RO from col above 
df5$RO_downf <- paste("RO=", df5$RO_downf, sep="")
df5$RO_downf <- paste(df5$RO_downf, ";", sep="")

#in new vcf - replace AO values here with AO from col above 
df5$AO_down <- paste("AO=", df5$AO_down, sep="")
df5$AO_down <- paste(df5$AO_down, ";", sep="")

#easier packacge to write to
library(vcfR)
library(stringr)

vcf <- read.vcfR("../475_sam_dupl_rm_no5_plasma6mo.vcf")
vcf <- read.vcfR("../475_sam_dupl_rm_no5_urine6mo.vcf")


#1 - cutdown to length of filter file above 
vcf@fix=vcf@fix[vcf@fix[,'POS'] %in% df5$start,]
#this works - but it does not change the rowID - no way to link @gt lines and @fix lines 
#apparently you dont need to drop the corresponding @gt lines just writes it out correct - cool 
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], "RO=(.*?);", as.character(df5$RO_downf))
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], "AO=(.*?);", as.character(df5$AO_down))



#will output it compressed regardless of filename
write.vcf(vcf, file="good_variants_B103plasma6mo.100.vcf.gz")


write.vcf(vcf, file="good_variants_B103urine6mo.100.vcf.gz")

#########1000x downsampling 
#so using df3 as the true number of visible sites - how many of them are polymorphic? 
#this is where we can then use qual etc 
df=df3[df3$DP > 100,]
df=df[df$QUAL > 0,]



#cut down to 100 frequency
df$RO_down <- round(1000*(df$RO/df$DP))
df$AO1_down <- round(1000*(as.numeric(df$AO1)/df$DP))
df$AO2_down <- round(1000*(as.numeric(df$AO2)/df$DP))
df$AO3_down <- round(1000*(as.numeric(df$AO3)/df$DP))

#get rid of sites where ref is greater than 98 - if 99 other sites cant be greater than 2%
df=df[df$RO_down <= 980,]
#get rid of sites where no alt is at least greater than 2% but also less than 99
df=df[df$AO1_down <=980 & df$AO1_down > 19 | df$AO2_down <=980 & df$AO2_down > 19 | df$AO3_down <=980 & df$AO3_down > 19,] 
df=df[complete.cases(df$Row.names),]

##THIS - THIS IS THE RIGHT IDEA 
#NEEDS TO BE AT LEAST TWO SITES MEETTHING THESE CONDITIONS - NOT OROROR 
#OTHERWISE ITS NOT TRULY POLYMORPHIC 
#SOMETHINGS BEING WEIRD WITH THE NA THO
df=as.data.frame(df, stringsAsFactors = FALSE)
#perfect! just throws a warning not an error!!
df[is.na(df)] <- 0

###whoops whoops cant od this 
#all those sites still need to be supported by at least one read one ither sites and at that site 
##tris=subset(df5,  rowSums(df5[65:68] >= 2) == 3)
#PREV CHUNK WAS FOR TOTAL THIS NOW BREAKS DOWN BY BI/TRI/QUAD ALLELEIC
# == INSTEAD OF >=
df5 <- df[rowSums(cbind((df$RO_down >= 20 & df$SRR > 0 & df$SRF > 0), (df$AO1_down >= 20 & df$SAR1 > 0 & df$SAF1 > 0), 
                        (df$AO2_down >= 20 & df$SAR2 > 0 & df$SAF2 > 0),  
                        (df$AO3_down >= 20 & df$SAR3 > 0 & df$SAF3 > 0))) == 2,] #CHANGE HERE FOR BIALLEICS (2) #CHANGE HERE FOR TRIALLEICS (3) #CHANGE HERE FOR QUADALLEICS (4) 
df5=df5[complete.cases(df5$Row.names),]
df5=subset(df5, ((RO_down >= 20) & rowSums(cbind((df5$AO1_down >= 20 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                                (df5$AO2_down >= 20 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                                (df5$AO3_down >= 20 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 1) | #CHANGE HERE FOR BIALLEICS (1) #CHANGE HERE FOR TRIALLEICS (2) #CHANGE HERE FOR QUADALLEICS (3) 
             ((RO_down < 20) & rowSums(cbind((df5$AO1_down >= 20 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                            (df5$AO2_down >= 20 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                            (df5$AO3_down >= 20 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 2)) #CHANGE HERE FOR BIALLEICS (2) #CHANGE HERE FOR TRIALLEICS (3) #CHANGE HERE FOR QUADALLEICS (4) 
df5=df5[complete.cases(df5$Row.names),]

#pick the largest of the AOs as the alt
#100-that number becomes RO
#should account for if there's multiple AOs as highest
library(dplyr)
df5=df5 %>% 
  mutate(AO_down = pmax(AO1_down, AO2_down, AO3_down))

df5$RO_downf=1000-df5$AO_down

#in new vcf - replace RO values here with RO from col above 
df5$RO_downf <- paste("RO=", df5$RO_downf, sep="")
df5$RO_downf <- paste(df5$RO_downf, ";", sep="")

#in new vcf - replace AO values here with AO from col above 
df5$AO_down <- paste("AO=", df5$AO_down, sep="")
df5$AO_down <- paste(df5$AO_down, ";", sep="")

#easier packacge to write to
library(vcfR)
library(stringr)

vcf <- read.vcfR("../475_sam_dupl_rm_no5_plasma6mo.vcf")
vcf <- read.vcfR("../475_sam_dupl_rm_no5_urine6mo.vcf")


#1 - cutdown to length of filter file above 
vcf@fix=vcf@fix[vcf@fix[,'POS'] %in% df5$start,]
#this works - but it does not change the rowID - no way to link @gt lines and @fix lines 
#apparently you dont need to drop the corresponding @gt lines just writes it out correct - cool 
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], "RO=(.*?);", as.character(df5$RO_downf))
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], "AO=(.*?);", as.character(df5$AO_down))

#will output it compressed regardless of filename
write.vcf(vcf, file="good_variants_B103plasma6mo.1000.vcf.gz")
write.vcf(vcf, file="good_variants_B103urine6mo.1000.vcf.gz")