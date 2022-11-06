library(gplots)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

#Read in coca datasets and remove unnecessary columns 
#(only need ID, hit, evalue, foldchange, length, RPKM)

folder_path <- paste0("/Users/Gabriel/Desktop/Prashanth/Labwork/",
"Erythroxylaceae project/Sequences/Batch 2 - February 2021/Complete_dataset/")

ecoca_48 <- read.csv(paste(folder_path, "E_coca_48.csv", sep=""),header=TRUE)
ecoca_48 <- ecoca_48[,c(1,2,3,6,18)]

ecoca_113 <- read.csv(paste(folder_path, "E_coca_113.csv", sep=""),header=TRUE)
ecoca_113 <- ecoca_113[,c(1,2,3,6,18)]

ecoca_124 <- read.csv(paste(folder_path, "E_coca_124.csv", sep=""),header=TRUE)
ecoca_124 <- ecoca_124[,c(1,2,3,6,18)]

ecoca_209 <- read.csv(paste(folder_path, "E_coca_209.csv", sep=""),header=TRUE)
ecoca_209 <- ecoca_209[,c(1,2,3,6,18)]

ecoca_228 <- read.csv(paste(folder_path, "E_coca_228.csv", sep=""),header=TRUE)
ecoca_228 <- ecoca_228[,c(1,2,3,6,18)]

ehond <- read.csv(paste(folder_path, "E_hond.csv", sep=""), header=TRUE)
ehond <- ehond[,c(1,2,3,6,18)]

#Apply evalue cutoff (<=1E-50): remove genes without strong BLAST matches

ecoca_48 <- ecoca_48[ecoca_48$evalue <= 1E-50,]
ecoca_113 <- ecoca_113[ecoca_113$evalue <= 1E-50,]
ecoca_124 <- ecoca_124[ecoca_124$evalue <= 1E-50,]
ecoca_209 <- ecoca_209[ecoca_209$evalue <= 1E-50,]
ecoca_228 <- ecoca_228[ecoca_228$evalue <= 1E-50,]
ehond <- ehond[ehond$evalue <= 1E-50,]

#Sort by RPKM

ecoca_48 <- ecoca_48[order(-ecoca_48$RPKM),]
ecoca_113 <- ecoca_113[order(-ecoca_113$RPKM),]
ecoca_124 <- ecoca_124[order(-ecoca_124$RPKM),]
ecoca_209 <- ecoca_209[order(-ecoca_209$RPKM),]
ecoca_228 <- ecoca_228[order(-ecoca_228$RPKM),]
ehond <- ehond[order(-ehond$RPKM),]

#Add numbering to duplicates so each hit name is unique

ecoca_48$hit <- make.unique(ecoca_48$hit)
ecoca_113$hit <- make.unique(ecoca_113$hit)
ecoca_124$hit <- make.unique(ecoca_124$hit)
ecoca_209$hit <- make.unique(ecoca_209$hit)
ecoca_228$hit <- make.unique(ecoca_228$hit)
ehond$hit <- make.unique(ehond$hit)

#Make master list of BLAST match (hit) names and remove duplicates

hits <- c(ecoca_48$hit, ecoca_113$hit, ecoca_124$hit, ecoca_209$hit, 
          ecoca_228$hit, ehond$hit)
hits <- hits[!duplicated(hits)]

#Construct matrix of RPKM values across E. coca and E. hond for each hit 
#(assign RPKM = 1 if gene not found)

profile <- data.frame(hits)

profile$ecoca_48 <- ecoca_48$RPKM[match(hits,ecoca_48$hit)]
profile$ecoca_113 <- ecoca_113$RPKM[match(hits,ecoca_113$hit)]
profile$ecoca_124 <- ecoca_124$RPKM[match(hits,ecoca_124$hit)]
profile$ecoca_209 <- ecoca_209$RPKM[match(hits,ecoca_209$hit)]
profile$ecoca_228 <- ecoca_228$RPKM[match(hits,ecoca_228$hit)]
profile$ehond <- ehond$RPKM[match(hits,ehond$hit)]

profile[is.na(profile)] <- 1
profile_data <- data.matrix(log2(profile[,-1]))
rownames(profile_data) <- profile$hits

#Generate hierarchical clustering / heatmap

pheatmap(profile_data, cluster_rows=TRUE, cluster_cols=FALSE, scale="row",
         fontsize_row = 3)

##################################################################################################################
#Normalize RPKM values to housekeeping gene: 20DD_coca-20110131_13971 (highest-expressed in nearly all samples)

ecoca_48$norm_RPKM <- ecoca_48$RPKM / ecoca_48$RPKM[match("20DD_coca-20110131_13971",ecoca_48$hit)]
ecoca_113$norm_RPKM <- ecoca_113$RPKM / ecoca_113$RPKM[match("20DD_coca-20110131_13971",ecoca_113$hit)]
ecoca_124$norm_RPKM <- ecoca_124$RPKM / ecoca_124$RPKM[match("20DD_coca-20110131_13971",ecoca_124$hit)]
ecoca_209$norm_RPKM <- ecoca_209$RPKM / ecoca_209$RPKM[match("20DD_coca-20110131_13971",ecoca_209$hit)]
ecoca_228$norm_RPKM <- ecoca_228$RPKM / ecoca_228$RPKM[match("20DD_coca-20110131_13971",ecoca_228$hit)]
ehond$norm_RPKM <- ehond$RPKM / ehond$RPKM[match("20DD_coca-20110131_13971",ehond$hit)]

