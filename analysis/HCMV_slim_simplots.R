#ofiicla official final final figure
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
#have to filter by congenital and plasma and 100 vs 1000X sampling 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/")

data=read.csv("100bicombos2_all.csv",header=FALSE)
data=read.csv("50bicombos_2_notF2.csv",header=FALSE)


data2 = data %>%
  extract(V1, into = c("First", "Second"), "^([^.]+)\\.(.*)")
data2=data2[rowSums(is.na(data2)) != ncol(data2),]

data2[c('combo', 'snps','filter')] <- str_split_fixed(data2$Second,  "\\s+", 3)
data2[c('compartment_subsample', 'compartment')] <- str_split_fixed(data2$First, 'congenital', 2)
data2$compartment_subsample=gsub("\\s+", "", data2$compartment_subsample)
data2$combo <- sub(".[^.]+$", "", data2$combo)
data2$combo <- sub(".[^.]+$", "", data2$combo)
data2$combo <- sub(".[^.]+$", "", data2$combo)

data3=data2[data2$filter == "good", ]
data3=data3[!grepl("tri", data3$First),]



#barchart only of frequency of combos passing 
#need to get one subset of combos 
data4=data3[data3$compartment_subsample == "plasma_100_bi", ]

ggplot() + 
  geom_bar(data=data4, aes(combo), position="dodge") +
  theme(axis.text.x = element_text(angle = 40, vjust=1, hjust=1, size=14),
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) + 
  
  scale_x_discrete(labels = str_wrap(c("DFE1, unscaled μ=variable, unscaled θ=variable, skew=0.07", 
                                       "DFE1, unscaled μ=variable, unscaled θ=9.8e-06, skew=0.07", 
                                       "DFE1, unscaled μ=variable, unscaled θ=9.8e-07, skew=0.07", 
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.07", 
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.1",
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-08, skew=0.07", 
                                       "DFE2, unscaled μ=variable, unscaled θ=variable, skew=0.01",  
                                       "DFE2, unscaled μ=variable, unscaled θ=9.8e-06, skew=0.01", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=variable, skew=0.01",  
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=varible, skew=0.07", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=variable, skew=0.1", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-06, skew=0.07", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-06, skew=0.1", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.01", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.07"
                                       
  ),
  width=10))

# grouped boxplot
ggplot(data3, aes(x=combo, y=as.numeric(snps), fill=compartment_subsample)) + 
  geom_boxplot() + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 90, vjust=1, hjust=1,size = 13), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  geom_hline(aes(yintercept=431), color="salmon") + 
  geom_hline(aes(yintercept=304), color="#7CAE00") + 
  geom_hline(aes(yintercept=431), color="#00BFC4") +
  geom_hline(aes(yintercept=335), color="#C77CFF") +
  
  ggtitle("Filtered on +/- 100 biallelic and +/ 10 triallelic snps") +
  scale_x_discrete(labels = str_wrap(c("DFE1, unscaled μ=variable, unscaled θ=variable, skew=0.07", 
                                       "DFE1, unscaled μ=variable, unscaled θ=9.8e-06, skew=0.07", 
                                       "DFE1, unscaled μ=variable, unscaled θ=9.8e-07, skew=0.07", 
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.07", 
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.1",
                                       "DFE1, unscaled μ=2.0e-06, unscaled θ=9.8e-08, skew=0.07", 
                                       "DFE2, unscaled μ=variable, unscaled θ=variable, skew=0.01",  
                                       "DFE2, unscaled μ=variable, unscaled θ=9.8e-06, skew=0.01", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=variable, skew=0.01",  
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=varible, skew=0.07", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=variable, skew=0.1", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-06, skew=0.07", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-06, skew=0.1", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.01", 
                                       "DFE2, unscaled μ=2.0e-06, unscaled θ=9.8e-07, skew=0.07"
                                       
  ),
  width=10))


################################################################################################################################################
#F2 simulations - formatted slightly differently bc of replicate number in filename so reading in seperately
#have to filter by congenital and plasma and 100 vs 1000X sampling 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/")

#Not using the backinbibin strategy
data=read.csv("combos_100bi_nobibin.csv",header=FALSE)

#using the same backinbibin strategy implemented in the NON F2 simulations
data=read.csv("combos_100bi_bibin.csv",header=FALSE)

#using backinbibin strategy and lowering thresh to 50 - F2 admix sims
data=read.csv("combos_50_bibin.csv",header=FALSE)

data2 = data %>%
  extract(V1, into = c("First", "Second"), "^([^.]+)\\.(.*)")
data2=data2[rowSums(is.na(data2)) != ncol(data2),]

data2[c('combo', 'snps','filter')] <- str_split_fixed(data2$Second,  "\\s+", 3)

data2$combo=sub(".100.output.ms", "", data2$combo)
data2$combo=sub("./F2_rescaled/congenital_plasma.", "", data2$combo)
data2$combo=sub(".[^.]+$", "", data2$combo)



data3=data2[data2$filter == "good", ]


#barchart only of frequency of combos passing 
#need to get one subset of combos 
data4=data3[data3$First == "  plasma_100 ", ]

ggplot() + 
  geom_bar(data=data4, aes(combo), position="dodge") +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, size=8),
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank())

group.colors=c("salmon", "#7CAE00", "#00BFC4", "#C77CFF")

group.colors=c("#00BFC4","salmon","#C77CFF","#7CAE00")


# grouped boxplot
ggplot(data3, aes(x=combo, y=as.numeric(snps), fill=First)) + 
  geom_boxplot() + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 60, vjust=1, hjust=1,size = 8), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  geom_hline(aes(yintercept=433), color="salmon") + 
  geom_hline(aes(yintercept=304), color="#7CAE00") + 
  geom_hline(aes(yintercept=431), color="#00BFC4") +
  geom_hline(aes(yintercept=335), color="#C77CFF") +
  scale_fill_manual(values=group.colors)

