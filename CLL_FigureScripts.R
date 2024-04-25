# clean-up
rm(list=ls())
# load libraries
library("BoolNet") #we used version 2.1.7
library("qgraph") #we used version 1.9.3
library("igraph") #we used version 1.3.5
library("pROC") #we used version 1.18.0
library("ggplot2") #we used version 3.4.1
library("ggpubr") #we used version 0.6.0
library("poweRlaw") #we used version 0.70.6
library("readr") # we use version 2.1.5
library("viridis") # we use version 0.6.4
library("tidyverse") # we used version 2.0.0
# set working directory
setwd("./")
#load required custom functions
source("Functions/Required_Functions_CLL.R")

#load network and grouping
CLL <- loadNetwork("CLL_Model.txt")
groupingCLL<-list(class=rev(c("", "","", "", "" )),
                  index= rev(list(rev(c("BCR", "CD5","SYK", "LYN", "ZAP70", "BTK")),
                                  rev(c("PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2")),   
                                  rev(c("SHP1", "PTEN", "NFAT", "Anergy")), 
                                  rev(c("p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis")),
                                  rev(c("CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase")))))

CLL_TME<- loadNetwork("CLL_Model_Microenvironment.txt")
groupingCLL_TME<-list(class=rev(c("", "","", "", "" )),
                  index= rev(list(rev(c("BCR", "CD5", "TME","SYK", "LYN", "ZAP70", "BTK")),
                                  rev(c("PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2")),   
                                  rev(c("SHP1", "PTEN", "NFAT", "Anergy")), 
                                  rev(c("p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis")),
                                  rev(c("CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase")))))


##########################################################################################################
####################################### Figure 2 #########################################################
##########################################################################################################

#Figure 2A-D

#Test for scale-freeness
if(!require("poweRlaw",character.only = TRUE)) install.packages("poweRlaw")
library("poweRlaw",character.only = TRUE)

readNetwork <- function(nwFileName, nwName) {
  inFile <- paste(nwFileName,sep="")
  return(read.csv(inFile, stringsAsFactors=FALSE))
}

simplify <- function(genes, rules){
  srules <- list()
  for (r  in 1:length(rules)) {
    rule <- gsub(" ", "", rules[r])
    rule <- gsub("&", ",",rule, fixed=TRUE)
    rule <- gsub("|", ",",rule, fixed=TRUE)
    rule <- gsub("(", "", rule, fixed=TRUE)
    rule <- gsub(")", "", rule, fixed=TRUE)
    rule <- gsub("!","", rule, fixed=TRUE)
    srules[r] <- list(which(genes %in% strsplit(rule,",")[[1]]))
  }
  return (srules)
}


countDegree <- function(network, name, noOfSims = 500, noOfThreads=8 ) {
  ccNw <- readNetwork(network, name)
  ccGenes <- ccNw[,1]
  ccRules <- ccNw[,2]
  sr <- simplify(ccGenes,ccRules)
  degs <- matrix(0,ncol=6, nrow=length(ccGenes), dimnames=list(ccGenes,list("out","inp","loop","total","Z (out)", "Z (total)")))
  degs[,1] <- sapply(1:length(ccGenes),function(g) { length(sr[[g]]) })
  degs[,2] <- c(tabulate(unlist(sr)), rep(0,length(ccGenes)-max(unlist(sr))))
  degs[,3] <- sapply(1:length(ccGenes), function(g) { if (g %in% sr [[g]]) {return(1)} else {return(0)} } )
  degs[,1] <- degs[,1] - degs[,3]
  degs[,2] <- degs[,2] - degs[,3]
  degs[,4] <-  degs[,1] + degs[,2] + degs[,3]
  degMean <- mean(degs[,1])
  degSd <- sd(degs[,1])
  degs[,5] <- round((degs[,1] - degMean)/degSd,2)
  degMean <- mean(degs[,4])
  degSd <- sd(degs[,4])
  degs[,6] <- round((degs[,4] - degMean)/degSd,2)
  degtable <- degs
  degtable[degtable==0] <- NA
  m <- degs[,4]
  m <- m[m>0]
  m_pl <- displ$new(m)
  est <- estimate_xmin(m_pl)
  if(is.na(est$xmin)) {
    print("failed to set model parameters")
    pvalue <- "failed to set model parameters"
  } else {
    m_pl$setXmin(est)
    bt_pl <- bootstrap_p(m_pl, no_of_sims=noOfSims, threads=noOfThreads)
    pvalue <- bt_pl$p
  }
  return(list(bootstrap=bt_pl, degrees=degs))
}
set.seed(10000)
network <- "CLL_Model_Dec23.txt"
Zscore <- countDegree(network,network, noOfSims=100, noOfThreads=4)
print(Zscore)

##scale freeness result
#$p, to be considered scale free your p>0.1
#[1] 0.47
#Note, that genes/proteins with a score Z(total)>2.5 are defined as hub nodes. 

########Figure on inteGraph With corresponding Sizes#######
conv2adjmat <- function(sbmlnet, inputcorrected = FALSE){
  "Converts an SBML object to adjacency matrix. 
  If inputcorrected = T, all input edges are set to zero. 
  A vertex v is said to be an input if it is only regulated by itself, 
  meaning the sum of column v in the adjacency matrix is one, with the only non-zero entry 
  being at position [v,v]."
  adjmat <- sapply(sbmlnet$interactions, function(gene) {v <- rep(0,length(sbmlnet$genes)); 
  v[gene$input] <- 1; return(v)})
  if (inputcorrected == TRUE){
    for (d in 1:dim(adjmat)[1]){
      if (adjmat[d,d] >= 1){adjmat[d,d] <- 1}
    }
    for (col in 1:dim(adjmat)[2]){
      if (sum(adjmat[,col]) == 1 & adjmat[col,col] == 1){
        adjmat[col,col] <-0
      }
    }
  }
  return(adjmat)
}

CLL_adjmat <- conv2adjmat(CLL)


###### Interaction graph based on Boolean functions with node size related to Z-scores
totaldegs_noDelays <- rep(NA, length(CLL$genes))
for (g in 1:length(CLL$genes)){
  totaldegs_noDelays[g] <- sum(CLL_adjmat[g,]) + sum(CLL_adjmat[,g])
}
hist(totaldegs_noDelays, breaks = seq(1,18,1), main = "Total degree distribution in network without time delays")

totaldegs_zscores <- rep(NA, length(CLL$genes))
for (g in 1:length(CLL$genes)){
  totaldegs_zscores[g] <- (totaldegs_noDelays[g] - mean(totaldegs_noDelays))/sd(totaldegs_noDelays)
}

g <- graph_from_adjacency_matrix(CLL_adjmat)

#the following lines are not required, the user can already upload the coordinates from tkplots availbla as .RDS
####need to make a layout for the graph --> tkplot(g)
#tkid<- tkplot(g, canvas.width = 1000, canvas.height = 1000, layout=l)
#the next comand is to be called when the layout is open
#l<- tkplot.getcoords(tkid)
#saveRDS(l, file="plotCLLLayout.RDS")


#####Coloring interaction graph
#proliferation Genes
greengenes <- c("CTNNB1", "cmyc", "p21","p27", "p16", "CyclinD", "CyclinE", "Rb", "E2F", "SPhase", "p15" ) 
#Anergy Genes
bluegenes <- c("CD5", "SHP1", "PTEN", "NFAT", "Anergy")
#BCR Signalosome
yellowgenes<- c("BCR", "LYN", "SYK", "ZAP70")
#Apoptosis Genes
redgenes<- c("MDM2","p53", "MCL1", "BCL2", "BIM", "p14", "Apoptosis")

#rest, i.e. delay nodes remain white
colorvec <- rep("#D6D6D6", length(CLL$genes))
colorvec[which(CLL$genes %in% redgenes)] <- "#E2A096"
colorvec[which(CLL$genes %in% bluegenes)] <- "#9BCDD6"
colorvec[which(CLL$genes %in% yellowgenes)] <- "#F7DD93"
colorvec[which(CLL$genes %in% greengenes)] <- "#B0DC00"

nodelabelsize <- (2*totaldegs_zscores+10)*0.05
l <- readRDS("./plotCLLLayout.RDS")
plot(g, layout=l, vertex.size=2*totaldegs_zscores+10, vertex.label.color="black",
     vertex.color=colorvec, edge.color="grey", vertex.label.cex=nodelabelsize,
     edge.arrow.size=0.3)

####Degrees histograms (first one is up above) 
#Total degree = sum up both row and col for each gene to get a vector          
totaldegs_noDelays <- rep(NA, length(CLL$genes))
for (g in 1:length(CLL$genes)){
  totaldegs_noDelays[g] <- sum(CLL_adjmat[g,]) + sum(CLL_adjmat[,g])
}
# hist(totaldegs_noDelays, breaks = seq(1,18,1), main = "Total degree distribution in network without time delays")

totaldegs_noDelays_zscores <- rep(NA, length(CLL$genes))
for (g in 1:length(CLL$genes)){
  totaldegs_noDelays_zscores[g] <- (totaldegs_noDelays[g] - mean(totaldegs_noDelays))/sd(totaldegs_noDelays)
} 

ggplot(as.data.frame(totaldegs_noDelays), aes(x=totaldegs_noDelays)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Total node degree distribution in network without time delays") +
  xlab("Total node degree") + ylab("Frequency") + 
  theme_bw() + scale_y_continuous(breaks=seq(0,10,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,18,1)) + 
  theme(panel.grid.minor = element_blank())

#In-degree = sum up cols across the adjmat to get a vector  
indegs_noDelays <- colSums(CLL_adjmat)
# hist(indegs_noDelays) 
ggplot(as.data.frame(indegs_noDelays), aes(x=indegs_noDelays)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Node in-degree distribution in network without time delays") +
  xlab("Node in-degree") + ylab("Frequency") + 
  theme_bw() + scale_y_continuous(breaks=seq(0,21,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,11,1)) + 
  theme(panel.grid.minor = element_blank())

#Out-degree = sum up rows across the adjmat to get a vector
outdegs_noDelays <- rowSums(CLL_adjmat)
# hist(outdegs_noDelays)

ggplot(as.data.frame(outdegs_noDelays), aes(x=outdegs_noDelays)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Node out-degree distribution in network without time delays") +
  xlab("Node out-degree") + ylab("Frequency") + 
  theme_bw() + scale_y_continuous(breaks=seq(0,19,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,13,1)) + 
  theme(panel.grid.minor = element_blank())

#totalDegree
ggplot(as.data.frame(totaldegs_noDelays), aes(x=totaldegs_noDelays)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Total degree distribution") +
  xlab("Node out-degree") + ylab("Frequency") + 
  theme_bw() + scale_y_continuous(breaks=seq(0,19,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,13,1)) + 
  theme(panel.grid.minor = element_blank())

###### Fig2F
###### - Stability assessment

set.seed(10000)
a<-loadNetwork(file = "CLL_Model_Dec23.txt")
Hamming<-modTestNetworkProperties(a, numRandomNets = 1000, 
                                  testFunction = "testTransitionRobustness",
                                  testFunctionParams = list(numSamples=1000),
                                  alternative="less")

#results for the Draft in testing stability of the network
#[1] "orig Result :"
#[1] 0.02012245
#[1] "random Result :"
#5% 
#0.0330602 

###### Fig2H
###### - CLL attractors
set.seed(100000)

groupingCLL<-list(class=rev(c("", "","", "", "" )),
                  index= rev(list(rev(c("BCR", "CD5","SYK", "LYN", "ZAP70", "BTK")),
                                  rev(c("PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2")),   
                                  rev(c("SHP1", "PTEN", "NFAT", "Anergy")), 
                                  rev(c("p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis")),
                                  rev(c("CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase")))))
#to fasten the simulations provide less startStates in the getattractors function, or use method = "sat.exhautive" instead 
attr_CLL <- getAttractors(CLL, startStates = 10000000)
plotAttractors(attr_CLL, allInOnePlot = T, grouping = groupingCLL, drawLegend = F ,reverse= T, title = "unperturbed CLL model", onColor = "#50A804", offColor = "steelblue4")


##########################################################################################################
####################################### Figure 3 #########################################################
##########################################################################################################

###### Fig3A
###### - Fused Attractor patterns CLL & RS

#loadmodels
CLL<- loadNetwork("CLL_Model.txt")
CLL_CDKN2a_TP53 <- fixGenes(CLL, c("p53","p16","p14","p15"), c(0,0,0,0))
CLL_AKT <- fixGenes(CLL, "AKT", 1)
CLL_NOTCH<- fixGenes(CLL, "NOTCH1",1)
CLL_NFAT<- fixGenes(CLL, "NFAT", 0)
CLL_TME<- loadNetwork("CLL_Model_Microenvironment.txt")
#simulatemodels
set.seed(10000)
#to fasten the simulations provide less startStates in the getattractors function, or use method = "sat.exhautive" instead 
a1 <- getAttractors(CLL, startStates = 10000000)
a2 <- getAttractors(CLL_CDKN2a_TP53, startStates = 10000000)
a3 <- getAttractors(CLL_AKT, startStates = 10000000)
a4 <- getAttractors(CLL_NOTCH, startStates = 10000000)
a5 <- getAttractors(CLL_NFAT, startStates = 10000000)
a10 <- getAttractors(CLL_TME, startStates = 10000000)

#getattractorsfused
CLL_Fused <-OverallActivities(a1)
CLL_CDKN2a_TP53_Fused <-OverallActivities(a2)
CLL_AKT_Fused <- as.data.frame(AttractorActivityWeightedBasin(a3, 1))
CLL_NOTCH_Fused <- as.data.frame(AttractorActivityWeightedBasin(a4, 1))
CLL_NFAT_Fused <- OverallActivities(a5)
CLL_TME_FUsed<- CLL_Fused <-OverallActivities(a10)

rownames(CLL_Fused)<-factor(rownames(CLL_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))

#plots
#CLL
a <- factor(rownames(CLL_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_Fused, aes(x = "rowSums(AllIn)", y = a, fill=CLL_Fused[,1])) +
  geom_tile() + ggtitle("CLL") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#CDKN2A_TP53
b <- factor(rownames(CLL_CDKN2a_TP53_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_CDKN2a_TP53_Fused, aes(x = "rowSums(AllIn)", y = b, fill=CLL_CDKN2a_TP53_Fused[,1])) +
  geom_tile() + ggtitle("CLL_CDKN2ab_TP53") +  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#AKT
c <- factor(rownames(CLL_AKT_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_AKT_Fused, aes(x = "rowSums(AllIn)", y = c, fill=CLL_AKT_Fused[,1])) +
  geom_tile() +ggtitle("CLL_AKT") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#NOTCH
d <- factor(rownames(CLL_NOTCH_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_NOTCH_Fused, aes(x = "rowSums(AllIn)", y = d, fill=CLL_NOTCH_Fused[,1])) +
  geom_tile() + ggtitle("CLL_NOTCH") +scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#NFAT
e <- factor(rownames(CLL_NFAT_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_NFAT_Fused, aes(x = "rowSums(AllIn)", y = e, fill=CLL_NFAT_Fused[,1])) +
  geom_tile() + ggtitle("CLL_Nfat") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))

l<- factor(rownames(CLL_TME_FUsed), levels = c("BCR", "CD5", "TME", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_TME_FUsed, aes(x = "rowSums(AllIn)", y = l, fill=CLL_TME_FUsed[,1])) +
  geom_tile() + ggtitle("CLL_TME") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))

###### Fig3B
###### - Attractor patterns CLL & RS

#load models and define conditions
#to fasten the simulations provide less startStates in the getattractors function, or use method = "sat.exhautive" instead 
CLL<- loadNetwork("CLL_Model.txt")
CLL_BCR_KO<- fixGenes(CLL, c("BCR"), c(0))
high_risk_CLL <- fixGenes(CLL, c("p53", "cmyc", "ZAP70"), c(0,1,1))
p15p14p15_ko <- fixGenes(CLL, c("p16","p14", "p15"), c(0,0,0))
p15p14p15_ko_TP53_KO <- fixGenes(CLL, c("p16","p53", "p14", "p15"), c(0,0,0,0))
p15p14p15_ko_TP53_KO_BCRKO <- fixGenes(CLL, c("p16","p53", "p14", "p15", "BCR"), c(0,0,0,0,0))
p15p14p15_ko_TP53_KO_MYC_KI <- fixGenes(CLL, c("p16","p53", "p14", "p15", "cmyc"), c(0,0,0,0,1))
AKT_KI <- fixGenes(CLL, c("AKT"), c(1)) 
NFAT_KO <- fixGenes(CLL, c("NFAT"), c(0)) 
CLL_TME<- loadNetwork("CLL_Model_Microenvironment.txt")

#simulate attractors and plot
set.seed(10000)
conditions <- list(CLL = CLL, CLL_BCR_KO = CLL_BCR_KO, high_risk_CLL = high_risk_CLL, p15p14p15_ko = p15p14p15_ko, p15p14p15_ko_TP53_KO = p15p14p15_ko_TP53_KO, p15p14p15_ko_TP53_KO_BCRKO = p15p14p15_ko_TP53_KO_BCRKO, p15p14p15_ko_TP53_KO_MYC_KI = p15p14p15_ko_TP53_KO_MYC_KI, AKT_KI= AKT_KI, NFAT_KO = NFAT_KO, CLL_TME = CLL_TME)

lapply(names(conditions), function(cond_name) {
  attr <- getAttractors(conditions[[cond_name]], startStates = 10000000)
  plotAttractors(attr, 
                 allInOnePlot = TRUE, 
                 grouping = groupingCLL, 
                 drawLegend = FALSE, 
                 reverse = TRUE, 
                 title = paste0(cond_name, " - CLL model"), 
                 onColor = "#50A804", 
                 offColor = "steelblue4")
})

####pie charts
#pdf("PieChartsFigure1B.pdf", width = 10, height =20)
#colorCoding
#316FAB --> anergy
#6EB944 --> proliferation
#A1CE82 --> cell cylce allert
#D66929 --> apoptosis
#unperturbed - CLL 
slices<- c(99.25, 0.07, 0.68)
lbls<- c("Anergy", "Proliferation", "Apoptosis") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#316FAB", "#6EB944", "#D66929"), main = "Unperturbed model")

#BCR KO 
slices<- c(100)
lbls<- c("apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#D66929"), main = "BCR KO")


#HIGH risk CLL KO 
slices<- c(16.7, 81.06, 2.24)
lbls<- c("Anergy", "cell cycle alert", "Proliferation") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#4E4EC2","#A1CE82", "#6EB944"), main = "HIGH risk CLL")
dev.off()

#p16/p14/p15
slices<- c(51.4, 0.73, 6.17,41.7)
lbls<- c("Anergy","Proliferation", "Apoptosis", "cell cycle alert") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#4E4EC2", "#6EB944", "#D66929","#A1CE82" ), main = "p16/p14/p15 KO")


#p16/p14/p15/p53 KO
slices<- c(99.69, 0.31)
lbls<- c("Proliferation", "Apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EB944", "#D66929" ), main = "p16/p14/p15/p53 KO")


#p16/p14/p15/p53/BCR KO
slices<- c(100)
lbls<- c( "Apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c( "#D66929" ), main = "p16/p14/p15/p53/BCR KO")

#AKT_NOTCH KI
slices<- c(100)
lbls<- c("Proliferation") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EB944"), main = "AKT_NOTCH KI")

#NFAT KO
slices<- c(99.37, 0.63)
lbls<- c("Proliferation", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EB944", "#D66929"), main = "NFAT KO")

#TME 
slices<- c(57.07, 0.37, 45.56)
lbls<- c("Anergy","Apoptosis","Proliferation") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#4E4EC2", "#D66929", "#6EB944"), main = "NOTCH KI")
#dev.off()

##########################################################################################################
####################################### Figure 4 #########################################################
##########################################################################################################

#single cell anaylses are provided as a separate Rmarkdown file

###### Fig4A
###### - CD5 levels on patient cohorts
CD5_cohort <- read_delim("CD5_Expression_Anonymous.csv",  delim = ",", escape_double = FALSE, trim_ws = TRUE)


#Binarization
rocCD5 <- roc(CD5_cohort$Entity, CD5_cohort$CD5, percent=TRUE,partial.auc.correct=FALSE,
              plot=TRUE,thresholds="best", 
              print.thres="best", print.auc = TRUE, main="CD5")
wilcox.test(CD5_cohort$CD5~CD5_cohort$Entity,alternative=c("two.sided"), paired=F, conf.level=0.95)
#pdf("CD5_Expression.pdf", width = 9, height =8.5)
ggplot(CD5_cohort, aes(x=Entity, y=CD5)) + 
  geom_boxplot(width = 0.3, fill="#009ACD") +
  geom_jitter(color="black", size=1.5, alpha=0.9, width = 0.03) +
  scale_fill_viridis(discrete = FALSE, alpha=0.06) +
  labs(x="",y="CD5 expression") +
  theme(
    legend.position="none",
    plot.title = element_text(size=20)) +
  theme(
    axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),  
    axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
    axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 1, face = "plain"),
    legend.position = "none")  + 
  geom_hline(yintercept = 230, linetype="dashed", color="red") +
  xlab("") + 
  ylab ("CD5 Expression")
#dev.off()

##########################################################################################################
####################################### Figure 5 #########################################################
##########################################################################################################

###### Fig 5c
###### - Input Output analysis

##required variables for all runs
CLL<- loadNetwork("CLL_Model.txt")

set.seed(10000)
inputActivities <- matrix(c(0,0,100,100), nrow = 2, ncol = 2)

#The following two lines are used to produce the input output analysis. Since the calculation is computationally intensive, we already provide the RData file contantining the simulation results.
#CLL_simresult1 <- InputOutputHeatmap(CLL_TME, inputs = c("NOTCH1", "TME"), plottedOutputGene = "BCR", inputActivities = inputActivities, nrStartStates = 250, nrTransitions = 500)
#save(CLL_simresult1, file = "CLL_IO-analysis_simResults.RData")

CLL_simresult1 <- get(load("CLL_IO-Notch_TP53.RData"))

#Script for producing the matrices for each gene
#pdf("CLL_IO-analysis_Microenv.pdf")
for (node in CLL$genes) {
  InputOutputHeatmap(CLL, inputs = c("NOTCH1", "p53"), plottedOutputGene = node, inputActivities = inputActivities, nrStartStates = 100, nrTransitions = 250,  simResultList=CLL_simresult1)
}
#dev.off()

###### Fig 5D
###### - Input Output analysis

CLL<- loadNetwork("CLL_Model.txt")
p141516p53 <- fixGenes(CLL, c("p16","p53", "p14", "p15"), c(0,0,0,0))
Aktki <- fixGenes(CLL, c("AKT"), c(1))
NFATko <- fixGenes(CLL, c("NFAT"), c(0))
#startstates for Cascades
CLLStartState2 <-generateState(CLL, specs= c(BCR=1), default = 0)
CLLStartState3 <-generateState(CLL, specs= c(BCR=1, AKT=1), default = 0)

#plot cascades
#pdf("CascadesToPhenotypes.pdf", width = 10, height =20)
CLLSeq1 <- plotSequence (CLL, startState = CLLStartState2, drawLegend = F, grouping = groupingCLL, title = "Unperturbed model", onColor = "#50A804", offColor = "steelblue4", reverse= T )
CLLSeq2 <- plotSequence (p141516p53, startState = CLLStartState2, drawLegend = F, grouping = groupingCLL, reverse= T, title = "CDKN2AB_p53_KO", onColor = "#50A804", offColor = "steelblue4" )
CLLSeq3 <- plotSequence (Aktki, startState = CLLStartState3, drawLegend = F, grouping = groupingCLL, reverse= T, title = "AKT_KI", onColor = "#50A804", offColor = "steelblue4" )
CLLSeq4 <- plotSequence (NFATko, startState = CLLStartState2, drawLegend = F, grouping = groupingCLL, reverse= T, title = "NFAT_KO", onColor = "#50A804", offColor = "steelblue4" )
#dev.off()









##########################################################################################################
####################################### Figure 6 #########################################################
##########################################################################################################


###### Fig 6A
###### - Fused attractors Tumor Drivers
#generate models
CLL<- loadNetwork("CLL_Model.txt")
CLL_BMI1_KI <- fixGenes(CLL, "BMI1", 1)
CLL_TP53_KO <- fixGenes(CLL, "p53", 0)
CLL_BMI1KI_TP53KO<- fixGenes(CLL, c("BMI1","p53"), c(1,0))
CLL_HighRisk <- fixGenes(CLL, c("p53", "cmyc","ZAP70"), c(0,1,1))
#simulatemodels
set.seed(10000)
a6 <- getAttractors(CLL_BMI1_KI, startStates = 10000000)
a7 <- getAttractors(CLL_TP53_KO, startStates = 10000000)
a8 <- getAttractors(CLL_BMI1KI_TP53KO, startStates = 10000000)
a9 <- getAttractors(CLL_HighRisk, startStates = 10000000)
#fuse attractors
CLL_BMI1KI_Fused <-OverallActivities(a6)
CLL_TP53KO_Fused <- OverallActivities(a7)
CLL_BMI1KITP53KO_Fused <- OverallActivities(a8)
CLL_HighRisk_Fused <- OverallActivities(a9)
#plots
#BMI1
f <- factor(rownames(CLL_BMI1KI_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_BMI1KI_Fused, aes(x = "rowSums(AllIn)", y = f, fill=CLL_BMI1KI_Fused[,1])) +
  geom_tile() + ggtitle("BMI1 KI fused") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#TP53
g <- factor(rownames(CLL_TP53KO_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_TP53KO_Fused, aes(x = "rowSums(AllIn)", y = g, fill=CLL_TP53KO_Fused[,1])) +
  geom_tile() + ggtitle("TP53 KO fused") +  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#BMI1 & TP53
h <- factor(rownames(CLL_BMI1KITP53KO_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_BMI1KITP53KO_Fused, aes(x = "rowSums(AllIn)", y = h, fill=CLL_BMI1KITP53KO_Fused[,1])) +
  geom_tile() + ggtitle("BMI1 and TP53 KO fused") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))
#HighRisk CLL
i <- factor(rownames(CLL_HighRisk_Fused), levels = c("BCR", "CD5", "SYK", "LYN", "ZAP70", "BTK", "PI3K", "AKT", "PP2A", "SET", "FOXO", "GSK3B", "PLCG2", "Ca", "DAG", "PKCdelta", "PKCbeta", "Ras", "ERK", "Ca_ERK", "p38MAPK", "NFkB", "Stim1", "NOTCH1", "S100a4", "BMI1", "CITED2", "SHP1", "PTEN", "NFAT", "Anergy", "p53","MDM2","MCL1",  "BIM", "BCL2", "Apoptosis", "CTNNB1", "cmyc", "p21", "p16", "p14", "p15", "p27", "CyclinD", "CyclinE", "Rb", "E2F", "Sphase"))
ggplot(CLL_HighRisk_Fused, aes(x = "rowSums(AllIn)", y = i, fill=CLL_HighRisk_Fused[,1])) +
  geom_tile() + ggtitle("High RIsk CLL fused") + scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = c("steelblue4",  "#4a916e", "#50A804"))


###### Fig 6C-D
###### - BMI1 expression + Ki67 quantification

#BMI_Evaluation <- read_delim("BMI1_Evaluation_Coort.csv", delim = ";", escape_double = FALSE, col_types = cols(BMI1 = col_number(), `TP53` = col_number(), Ki67 = col_number()), trim_ws = TRUE)
BMI_Evaluation <- read_delim("BMI1_Expression_Anonymous.csv", delim = ",", escape_double = FALSE, col_types = cols(BMI1 = col_number(), `TP53` = col_number(), Ki67 = col_number()), trim_ws = TRUE)

BMI_Evaluation <- BMI_Evaluation %>%
  group_by(TP53, Diagnose) %>%
  mutate(mean=mean(BMI1, na.rm = TRUE)) %>% 
  print

BMI_Evaluation <- BMI_Evaluation %>%
  mutate(Entity=case_when(Diagnose =="CLL" & TP53 == 0 ~ "CLL TP53 \u394", 
                          Diagnose =="CLL" & TP53 == 1 ~ "CLL TP53 WT", 
                          Diagnose =="RS" & TP53 == 0 ~ "RS TP53 \u394", 
                          Diagnose =="RS" & TP53 == 1 ~ "RS TP53 WT",
                          Diagnose =="akz. CLL" & TP53 == 1 ~ "akz_TP53WT",
                          Diagnose =="akz. CLL" & TP53 == 0 ~ "akz_TP53Mut")) %>% 
  print

#pdf("BMI1_Expression.pdf", width = 12, height =6)
ggplot(BMI_Evaluation, aes(x = Entity, y=BMI1)) +
  geom_boxplot(width = 0.3, fill="#009ACD") +
  scale_fill_viridis(discrete = FALSE, alpha=0.6) +
  geom_jitter(color="black", size=1.5, alpha=0.9, width = 0.02) +
  theme(
    legend.position="none",
    plot.title = element_text(size=20)) +
  theme(
    axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),  
    axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
    axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 1, face = "plain")) +
  xlab("") + 
  ylab ("BMI1 Expression")
#dev.off()

#pdf("Ki67_Expression.pdf", width = 12, height =6)
ggplot(BMI_Evaluation, aes(x = Entity, y=Ki67)) +
  geom_boxplot(width = 0.3, fill="#FFA700") +
  scale_fill_viridis(discrete = FALSE, alpha=0.6) +
  geom_jitter(color="black", size=1.5, alpha=0.9, width = 0.02) +
  theme(
    legend.position="none",
    plot.title = element_text(size=20)) +
  theme(
    axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),  
    axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
    axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 1, face = "plain")) +
  xlab("") + 
  ylab ("Ki67 in %")
#dev.off()
subset_CLLandRS <- subset(BMI_Evaluation, Diagnose != "akz. CLL")
wilcox.test(subset_CLLandRS$BMI1[subset_CLLandRS$Entity=="CLL TP53 \u394"], 
            subset_CLLandRS$BMI1[subset_CLLandRS$Entity=="RS TP53 \u394"],
            alternative = "less", exact=FALSE)
wilcox.test(subset_CLLandRS$BMI1[subset_CLLandRS$Entity=="CLL TP53 WT"], 
            subset_CLLandRS$BMI1[subset_CLLandRS$Entity=="RS TP53 WT"],
            alternative = "less", exact=FALSE)





###### Fig 6G
###### - DrugTarget Screening

CLL<- loadNetwork("CLL_Model.txt")
###Drug Screening###
p141516_p53_SET <- fixGenes(CLL, c("p16","p53", "p14", "p15", "SET"), c(0,0,0,0,0))
p141516_p53_PKCB <- fixGenes(CLL, c("p16","p53", "p14", "p15", "PKCbeta"), c(0,0,0,0,0))
p141516_p53_PKCB_SET <- fixGenes(CLL,c("p16","p53", "p14", "p15", "PKCbeta", "SET"), c(0,0,0,0,0,0))
AKT_SET <- fixGenes(CLL, c("AKT", "SET"), c(1,0))
AKT_PKCB <- fixGenes(CLL, c("AKT", "PKCbeta"), c(1,0))
AKT_PCKB_SET <- fixGenes(CLL, c("AKT", "PKCbeta", "SET"), c(1,0,0))
NFAT_SET <- fixGenes(CLL, c("NFAT", "SET"), c(0,0))
NFAT_PKCB <- fixGenes(CLL, c("NFAT", "PKCbeta"), c(0,0))
NFAT_PKCB_SET <- fixGenes(CLL, c("NFAT", "PKCbeta", "SET"), c(0,0,0))


set.seed(10000)
conditions2 <- list(p141516_p53_SET=p141516_p53_SET, p141516_p53_PKCB=p141516_p53_PKCB, p141516_p53_PKCB_SET=p141516_p53_PKCB_SET, AKT_SET=AKT_SET, AKT_PKCB=AKT_PKCB, AKT_PCKB_SET=AKT_PCKB_SET, NFAT_SET=NFAT_SET, NFAT_PKCB=NFAT_PKCB, NFAT_PKCB_SET=NFAT_PKCB_SET  )

lapply(names(conditions2), function(cond_name) {
  attr <- getAttractors(conditions2[[cond_name]], startStates = 10000000)
  plotAttractors(attr, 
                 allInOnePlot = TRUE, 
                 grouping = groupingCLL, 
                 drawLegend = FALSE, 
                 reverse = TRUE, 
                 title = paste0(cond_name, " - CLL model"), 
                 onColor = "#50A804", 
                 offColor = "steelblue4")
})


#piecharts
#p141516_p53_SET
slices<- c(98.37, 1.63)
lbls<- c("Proliferation", "Apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EB944", "#D66929"), main = "p141516_p53_SET")

#p141516_p53_PKCB
slices<- c(98.21, 0.02, 1.75)
lbls<- c("cell cycle alert", "Proliferation", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EF100","#6EB944", "#D66929" ), main = "p141516_p53_PKCB")

#p141516_p53_PKCB_SET
slices<- c(97.86, 2.14)
lbls<- c("quiescence", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#316FFB", "#D66929" ), main = "p141516_p53_PKCB_SET")

#AKT_SET
slices<- c(98.5, 1.5)
lbls<- c("quiescence", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#316FFB", "#D66929" ), main = "AKT_SET")

#AKT_PKCB
slices<- c(98.5, 1.5)
lbls<- c("quiescence", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#316FFB", "#D66929" ), main = "AKT_PKCB")

#AKT_PKCB_SET
slices<- c(98.16, 1.84)
lbls<- c("quiescence", "apoptosis") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#316FFB", "#D66929" ), main = "AKT_PKCB_SET")


#NFAT_SET/NFAT_PKCB/NFAT_SET_PKCB
slices<- c(100)
lbls<- c("Proliferation") #Cell cycle alert
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#6EB944" ), main = "NFAT_SET/NFAT_PKCB/NFAT_SET_PKCB")










