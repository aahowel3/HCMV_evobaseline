library(dplyr)
args = commandArgs(trailingOnly=TRUE)
vcf <- read.csv(args[1], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=try ( vcftotal %>%
  select(where(function(x) any(x %in% subsam))), silent=TRUE)


if(inherits(vcftotal2, "try-error"))
{
  con_plasma_100_bi=0
  con_plasma_100_tri=0
}
  if(class(vcftotal2) == "data.frame") 
  {
    
  library(janitor)
  vcftotal2=vcftotal2 %>%
  row_to_names(row_number = 1)

  bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]
v=colnames(vcftotal2)
 tab=table(v)
 vbad=names(tab[tab > 2])
 vvbad=bad[which(!bad %in% vbad)]
#this is not the way to do it - if there are no duplicates returns an empty bionly column
  trionly=vcftotal2[ , which(names(vcftotal2) %in% vvbad), drop=FALSE]
  bionly <- vcftotal2[, !colnames(vcftotal2) %in% bad, drop=FALSE]

  bionly=as.data.frame(lapply(bionly, as.numeric))
  bionly2percent = bionly[,colSums(bionly[, -1, drop=FALSE]) >= 2, drop=FALSE]
  bionly2percent = bionly2percent[,colSums(bionly2percent[, -1, drop=FALSE]) <= 98, drop=FALSE]

#for trialleics need to add duplicate columns together
#and then EACH MAF has to be over 2% to pass 
#but need to work around the mutation stacking issue in a single genome in a single position 
  trionly=as.data.frame(lapply(trionly, as.numeric))

  combinedtri= try(as.data.frame(mapply(paste0, trionly[, seq(1, length(trionly), 2)], trionly[, seq(2,length(trionly), 2)])) , silent = TRUE ) 
  if(inherits(combinedtri, "try-error"))
  {
    con_plasma_100_tri=0
    con_plasma_100_bi=length(bionly2percent)
    
  }
  if(class(combinedtri) == "data.frame") {
    #drop=FALSE is important if something is reduced down to a single column keep it as a datafrma and keep col name
    combinedtri[combinedtri == "11"] <- "01"
    truetri=combinedtri[, sapply(combinedtri, function(col) ( (sum(col =="01") >= 2) & (sum(col =="10") >= 2) & (sum(col =="00") >= 2) 
                                                              & (sum(col =="01")) <= 98) & (sum(col =="10") <= 98) & (sum(col =="00") <= 98) ), drop=FALSE] 
    combinedtriless <- combinedtri[, !colnames(combinedtri) %in% colnames(truetri), drop=FALSE]
    backinbibin=try( combinedtriless[, sapply(combinedtriless, function(col) sum((sum(col =="01") >= 2) | ( (sum(col =="10") >= 2) | (sum(col =="00") >= 2)) >= 2 ) & 
                                           sum((sum(col =="01") <= 98) | ( (sum(col =="10") <= 98) | (sum(col =="00") <= 98)) >= 2 )), drop=FALSE] , silent=TRUE)
    if(inherits( backinbibin, "try-error"))
    {
    con_plasma_100_bi=length(bionly2percent) 
    con_plasma_100_tri=length(truetri)
    }
    if(class(backinbibin) == "data.frame"){
    con_plasma_100_bi=length(bionly2percent) + length(backinbibin)
    con_plasma_100_tri=length(truetri)
    }
    
  }
  
  }


###############################################################################
vcf <- read.csv(args[2], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=try ( vcftotal %>%
                  select(where(function(x) any(x %in% subsam))), silent=TRUE)


if(inherits(vcftotal2, "try-error"))
{
  con_plasma_1000_bi=0
  con_plasma_1000_tri=0
}
if(class(vcftotal2) == "data.frame") 
{
  
  library(janitor)
  vcftotal2=vcftotal2 %>%
    row_to_names(row_number = 1)
  
  bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]
  v=colnames(vcftotal2)
 tab=table(v)
 vbad=names(tab[tab > 2])
 vvbad=bad[which(!bad %in% vbad)]
  #this is not the way to do it - if there are no duplicates returns an empty bionly column
  trionly=vcftotal2[ , which(names(vcftotal2) %in% vvbad), drop=FALSE]
  bionly <- vcftotal2[, !colnames(vcftotal2) %in% bad, drop=FALSE]
  
  bionly=as.data.frame(lapply(bionly, as.numeric))
  bionly2percent = bionly[,colSums(bionly[, -1, drop=FALSE]) >= 20, drop=FALSE]
  bionly2percent = bionly2percent[,colSums(bionly2percent[, -1, drop=FALSE]) <= 980, drop=FALSE]
  
  
  #for trialleics need to add duplicate columns together
  #and then EACH MAF has to be over 2% to pass 
  #but need to work around the mutation stacking issue in a single genome in a single position 
  trionly=as.data.frame(lapply(trionly, as.numeric))
  
  combinedtri= try(as.data.frame(mapply(paste0, trionly[, seq(1, length(trionly), 2)], trionly[, seq(2,length(trionly), 2)])) , silent = TRUE ) 
  if(inherits(combinedtri, "try-error"))
  {
    con_plasma_1000_tri=0
    con_plasma_1000_bi=length(bionly2percent)
    
  }
  if(class(combinedtri) == "data.frame") {
    #drop=FALSE is important if something is reduced down to a single column keep it as a datafrma and keep col name
    combinedtri[combinedtri == "11"] <- "01"
    truetri=combinedtri[, sapply(combinedtri, function(col) ( (sum(col =="01") >= 20) & (sum(col =="10") >= 20) & (sum(col =="00") >= 20) 
                                                              & (sum(col =="01")) <= 980) & (sum(col =="10") <= 980) & (sum(col =="00") <= 980) ), drop=FALSE] 
    combinedtriless <- combinedtri[, !colnames(combinedtri) %in% colnames(truetri), drop=FALSE]
    backinbibin=try ( combinedtriless[, sapply(combinedtriless, function(col) sum((sum(col =="01") >= 20) | ( (sum(col =="10") >= 20) | (sum(col =="00") >= 20)) >= 2 ) & 
                                           sum((sum(col =="01") <= 980) | ( (sum(col =="10") <= 980) | (sum(col =="00") <= 980)) >= 2 )), drop=FALSE] , silent=TRUE)
    if(inherits( backinbibin, "try-error"))
    {
      con_plasma_1000_bi=length(bionly2percent) 
      con_plasma_1000_tri=length(truetri)
    }
    if(class(backinbibin) == "data.frame"){
      con_plasma_1000_bi=length(bionly2percent) + length(backinbibin)
      con_plasma_1000_tri=length(truetri)
    }
    
  }
  
}
#############################################################################
vcf <- read.csv(args[3], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=try ( vcftotal %>%
                  select(where(function(x) any(x %in% subsam))), silent=TRUE)


if(inherits(vcftotal2, "try-error"))
{
  con_urine_100_bi=0
  con_urine_100_tri=0
}
if(class(vcftotal2) == "data.frame") 
{
  
  library(janitor)
  vcftotal2=vcftotal2 %>%
    row_to_names(row_number = 1)
  
  bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]
v=colnames(vcftotal2)
 tab=table(v)
 vbad=names(tab[tab > 2])
 vvbad=bad[which(!bad %in% vbad)]  
  #this is not the way to do it - if there are no duplicates returns an empty bionly column
  trionly=vcftotal2[ , which(names(vcftotal2) %in% vvbad), drop=FALSE]
  bionly <- vcftotal2[, !colnames(vcftotal2) %in% bad, drop=FALSE]
  
  bionly=as.data.frame(lapply(bionly, as.numeric))
  bionly2percent = bionly[,colSums(bionly[, -1, drop=FALSE]) >= 2, drop=FALSE]
  bionly2percent = bionly2percent[,colSums(bionly2percent[, -1, drop=FALSE]) <= 98, drop=FALSE]
  
  #for trialleics need to add duplicate columns together
  #and then EACH MAF has to be over 2% to pass 
  #but need to work around the mutation stacking issue in a single genome in a single position 
  trionly=as.data.frame(lapply(trionly, as.numeric))
  
  combinedtri= try(as.data.frame(mapply(paste0, trionly[, seq(1, length(trionly), 2)], trionly[, seq(2,length(trionly), 2)])) , silent = TRUE ) 
  if(inherits(combinedtri, "try-error"))
  {
    con_urine_100_tri=0
    con_urine_100_bi=length(bionly2percent)
    
  }
  if(class(combinedtri) == "data.frame") {
    #drop=FALSE is important if something is reduced down to a single column keep it as a datafrma and keep col name
    combinedtri[combinedtri == "11"] <- "01"
    truetri=combinedtri[, sapply(combinedtri, function(col) ( (sum(col =="01") >= 2) & (sum(col =="10") >= 2) & (sum(col =="00") >= 2) 
                                                              & (sum(col =="01")) <= 98) & (sum(col =="10") <= 98) & (sum(col =="00") <= 98) ), drop=FALSE] 
    combinedtriless <- combinedtri[, !colnames(combinedtri) %in% colnames(truetri), drop=FALSE]
    backinbibin=try ( combinedtriless[, sapply(combinedtriless, function(col) sum((sum(col =="01") >= 2) | ( (sum(col =="10") >= 2) | (sum(col =="00") >= 2)) >= 2 ) & 
                                           sum((sum(col =="01") <= 98) | ( (sum(col =="10") <= 98) | (sum(col =="00") <= 98)) >= 2 ) ), drop=FALSE], silent = TRUE)
    if(inherits( backinbibin, "try-error"))
    {
      con_urine_100_bi=length(bionly2percent) 
      con_urine_100_tri=length(truetri)
    }
    if(class(backinbibin) == "data.frame"){
      con_urine_100_bi=length(bionly2percent) + length(backinbibin)
      con_urine_100_tri=length(truetri)
    }
    
  }
  
}
###################################################################
vcf <- read.csv(args[4], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=try ( vcftotal %>%
                  select(where(function(x) any(x %in% subsam))), silent=TRUE)


if(inherits(vcftotal2, "try-error"))
{
  con_urine_1000_bi=0
  con_urine_1000_tri=0
}
if(class(vcftotal2) == "data.frame") 
{
  
  library(janitor)
  vcftotal2=vcftotal2 %>%
    row_to_names(row_number = 1)
  
  bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]
  v=colnames(vcftotal2)
 tab=table(v)
 vbad=names(tab[tab > 2])
 vvbad=bad[which(!bad %in% vbad)]
  #this is not the way to do it - if there are no duplicates returns an empty bionly column
  trionly=vcftotal2[ , which(names(vcftotal2) %in% vvbad), drop=FALSE]
  bionly <- vcftotal2[, !colnames(vcftotal2) %in% bad, drop=FALSE]
  
  bionly=as.data.frame(lapply(bionly, as.numeric))
  bionly2percent = bionly[,colSums(bionly[, -1, drop=FALSE]) >= 20, drop=FALSE]
  bionly2percent = bionly2percent[,colSums(bionly2percent[, -1, drop=FALSE]) <= 980, drop=FALSE]
  
  #for trialleics need to add duplicate columns together
  #and then EACH MAF has to be over 2% to pass 
  #but need to work around the mutation stacking issue in a single genome in a single position 
  trionly=as.data.frame(lapply(trionly, as.numeric))
  
  combinedtri= try(as.data.frame(mapply(paste0, trionly[, seq(1, length(trionly), 2)], trionly[, seq(2,length(trionly), 2)])) , silent = TRUE ) 
  if(inherits(combinedtri, "try-error"))
  {
    con_urine_1000_tri=0
    con_urine_1000_bi=length(bionly2percent)
    
  }
  if(class(combinedtri) == "data.frame") {
    #drop=FALSE is important if something is reduced down to a single column keep it as a datafrma and keep col name
    combinedtri[combinedtri == "11"] <- "01"
    truetri=combinedtri[, sapply(combinedtri, function(col) ( (sum(col =="01") >= 20) & (sum(col =="10") >= 20) & (sum(col =="00") >= 20) 
                                                              & (sum(col =="01")) <= 980) & (sum(col =="10") <= 980) & (sum(col =="00") <= 980) ), drop=FALSE] 
    combinedtriless <- combinedtri[, !colnames(combinedtri) %in% colnames(truetri), drop=FALSE]
    backinbibin=try( combinedtriless[, sapply(combinedtriless, function(col) sum((sum(col =="01") >= 20) | ( (sum(col =="10") >= 20) | (sum(col =="00") >= 20)) >= 2 ) & 
                                           sum((sum(col =="01") <= 980) | ( (sum(col =="10") <= 980) | (sum(col =="00") <= 980)) >= 2 )), drop=FALSE] , silent=TRUE)
    if(inherits( backinbibin, "try-error"))
    {
      con_urine_1000_bi=length(bionly2percent) 
      con_urine_1000_tri=length(truetri)
    }
    if(class(backinbibin) == "data.frame"){
      con_urine_1000_bi=length(bionly2percent) + length(backinbibin)
      con_urine_1000_tri=length(truetri)
    }
    
  }
  
}



df = data.frame(
  column2 = c("plasma_100_bi", "plasma_1000_bi", "plasma_100_tri","plasma_1000_tri",
              "urine_100_bi", "urine_1000_bi", "urine_100_tri","urine_1000_tri"),
  column3 = c(args[1],args[1],args[1],args[1],args[1],args[1],args[1],args[1]),
  column4 = c(con_plasma_100_bi, con_plasma_1000_bi, con_plasma_100_tri,con_plasma_1000_tri,
              con_urine_100_bi, con_urine_1000_bi, con_urine_100_tri,con_urine_1000_tri))

df$column1=ifelse(((431 - 100) <= con_plasma_100_bi &  con_plasma_100_bi <= (431 + 100)) &
                    ((304 - 100) <= con_plasma_1000_bi & con_plasma_1000_bi <= (304 + 100)) &
                    ((431 - 100) <= con_urine_100_bi & con_urine_100_bi <= (431 + 100)) &
                    ((335 - 100) <= con_urine_1000_bi & con_urine_1000_bi <= (335 + 100)), "good", "bad")

###if you want it to be evaluted on tris and bis
df$column1=ifelse(((431 - 10) <= con_plasma_100_bi &  con_plasma_100_bi <= (431 + 10)) &
                    ((304 - 10) <= con_plasma_1000_bi & con_plasma_1000_bi <= (304 + 10)) &
                    ((431 - 10) <= con_urine_100_bi & con_urine_100_bi <= (431 + 10)) &
                    ((335 - 10) <= con_urine_1000_bi & con_urine_1000_bi <= (335 + 10)) &
                    ((9 - 10) <= con_plasma_100_tri &  con_plasma_100_tri <= (9 + 10)) &
                    ((1 - 10) <= con_plasma_1000_tri & con_plasma_1000_tri <= (1 + 10)) &
                    ((8 - 10) <= con_urine_100_tri & con_urine_100_tri <= (8 + 10)) &
                    ((6 - 10) <= con_urine_1000_tri & con_urine_1000_tri <= (6 + 10)), "good", "bad")

names(df) <- NULL
options(width=160)
print(df,row.names=FALSE)
