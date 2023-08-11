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

    colnames(bionly2percent)<-  sub("X", "", colnames(bionly2percent))
    x1=t(as.data.frame(colnames(bionly2percent)))
    write( "//", file = args[2], append = T)
    cat( "segsites:", length(bionly2percent), file = args[2], append = T)
    x1[1]=paste("\npositions:", x1[1])
    write.table( x1,
                 file = args[2], 
                 append = T,
                 sep = " ",
                 row.names = F,
                 
                 col.names = F,
                 na="",
                 quote = F)
    write.table( bionly2percent,
                 file = args[2], 
                 append = T,
                 sep = "",
                 row.names = F,
                 col.names = F,
                 na="",
                 quote = F)
    
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
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      bitri=cbind(bionly2percent, truetri)
collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)      
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[2], append = T)
      cat( "segsites:", length(bitri), file = args[2], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[2], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[2], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
      
    }
    if(class(backinbibin) == "data.frame"){
      con_plasma_100_bi=length(bionly2percent) + length(backinbibin)
      con_plasma_100_tri=length(truetri)
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      backinbibin[backinbibin == "00"] <- "0"
      backinbibin[backinbibin == "01"] <- "1"
      backinbibin[backinbibin == "10"] <- "1"
      
      bitri=cbind(bionly2percent, truetri, backinbibin)
collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)      
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[2], append = T)
      cat( "segsites:", length(bitri), file = args[2], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[2], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[2], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
    }
    
  }
  

#######################################################################################################################
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
    
    colnames(bionly2percent)<-  sub("X", "", colnames(bionly2percent))
    x1=t(as.data.frame(colnames(bionly2percent)))
    write( "//", file = args[4], append = T)
    cat( "segsites:", length(bionly2percent), file = args[4], append = T)
    x1[1]=paste("\npositions:", x1[1])
    write.table( x1,
                 file = args[4], 
                 append = T,
                 sep = " ",
                 row.names = F,
                 
                 col.names = F,
                 na="",
                 quote = F)
    write.table( bionly2percent,
                 file = args[4], 
                 append = T,
                 sep = "",
                 row.names = F,
                 col.names = F,
                 na="",
                 quote = F)
    
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
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      bitri=cbind(bionly2percent, truetri)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[4], append = T)
      cat( "segsites:", length(bitri), file = args[4], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[4], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[4], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
      
    }
    if(class(backinbibin) == "data.frame"){
      con_plasma_100_bi=length(bionly2percent) + length(backinbibin)
      con_plasma_100_tri=length(truetri)
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      backinbibin[backinbibin == "00"] <- "0"
      backinbibin[backinbibin == "01"] <- "1"
      backinbibin[backinbibin == "10"] <- "1"
      
      bitri=cbind(bionly2percent, truetri, backinbibin)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[4], append = T)
      cat( "segsites:", length(bitri), file = args[4], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[4], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[4], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
    }
    
  }
  
  
###############################################################################################################################
  vcf <- read.csv(args[5], header=FALSE)
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
    
    colnames(bionly2percent)<-  sub("X", "", colnames(bionly2percent))
    x1=t(as.data.frame(colnames(bionly2percent)))
    write( "//", file = args[6], append = T)
    cat( "segsites:", length(bionly2percent), file = args[6], append = T)
    x1[1]=paste("\npositions:", x1[1])
    write.table( x1,
                 file = args[6], 
                 append = T,
                 sep = " ",
                 row.names = F,
                 
                 col.names = F,
                 na="",
                 quote = F)
    write.table( bionly2percent,
                 file = args[6], 
                 append = T,
                 sep = "",
                 row.names = F,
                 col.names = F,
                 na="",
                 quote = F)
    
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
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      bitri=cbind(bionly2percent, truetri)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[6], append = T)
      cat( "segsites:", length(bitri), file = args[6], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[6], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[6], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
      
    }
    if(class(backinbibin) == "data.frame"){
      con_plasma_100_bi=length(bionly2percent) + length(backinbibin)
      con_plasma_100_tri=length(truetri)
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      backinbibin[backinbibin == "00"] <- "0"
      backinbibin[backinbibin == "01"] <- "1"
      backinbibin[backinbibin == "10"] <- "1"
      
      bitri=cbind(bionly2percent, truetri, backinbibin)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[6], append = T)
      cat( "segsites:", length(bitri), file = args[6], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[6], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[6], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
    }
    
  }
  
############################################################################################################
  vcf <- read.csv(args[7], header=FALSE)
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
    
    colnames(bionly2percent)<-  sub("X", "", colnames(bionly2percent))
    x1=t(as.data.frame(colnames(bionly2percent)))
    write( "//", file = args[8], append = T)
    cat( "segsites:", length(bionly2percent), file = args[8], append = T)
    x1[1]=paste("\npositions:", x1[1])
    write.table( x1,
                 file = args[8], 
                 append = T,
                 sep = " ",
                 row.names = F,
                 
                 col.names = F,
                 na="",
                 quote = F)
    write.table( bionly2percent,
                 file = args[8], 
                 append = T,
                 sep = "",
                 row.names = F,
                 col.names = F,
                 na="",
                 quote = F)
    
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
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      bitri=cbind(bionly2percent, truetri)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[8], append = T)
      cat( "segsites:", length(bitri), file = args[8], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[8], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[8], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
      
    }
    if(class(backinbibin) == "data.frame"){
      con_plasma_100_bi=length(bionly2percent) + length(backinbibin)
      con_plasma_100_tri=length(truetri)
      
      truetri[truetri == "00"] <- "0"
      truetri[truetri == "01"] <- "1"
      truetri[truetri == "10"] <- "1"
      backinbibin[backinbibin == "00"] <- "0"
      backinbibin[backinbibin == "01"] <- "1"
      backinbibin[backinbibin == "10"] <- "1"
      
      bitri=cbind(bionly2percent, truetri, backinbibin)
      collate <- function(x){
        x %>% select(sort(colnames(bitri))) 
      }
      bitri=collate(bitri)
      colnames(bitri)<-  sub("X", "", colnames(bitri))
      x1=t(as.data.frame(colnames(bitri)))
      write( "//", file = args[8], append = T)
      cat( "segsites:", length(bitri), file = args[8], append = T)
      x1[1]=paste("\npositions:", x1[1])
      write.table( x1,
                   file = args[8], 
                   append = T,
                   sep = " ",
                   row.names = F,
                   
                   col.names = F,
                   na="",
                   quote = F)
      write.table( bitri,
                   file = args[8], 
                   append = T,
                   sep = "",
                   row.names = F,
                   col.names = F,
                   na="",
                   quote = F)
    }
    
  }
  
  
  



