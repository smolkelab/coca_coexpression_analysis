library(pheatmap)

#Read in coca datasets and remove unnecessary columns (only need ID, hit, evalue, foldchange, length, RPKM)

folder_path <- "/Users/Gabriel/Desktop/Prashanth/Labwork/Erythroxylaceae project/Sequences/Batch 2 - February 2021/Complete_dataset/"

ecoca_48 <- read.csv(paste(folder_path, "E_coca_48.csv", sep=""), header=TRUE)
ecoca_48 <- ecoca_48[,c(1,2,3,6,18)]

ecoca_113 <- read.csv(paste(folder_path, "E_coca_113.csv", sep=""), header=TRUE)
ecoca_113 <- ecoca_113[,c(1,2,3,6,18)]

ecoca_124 <- read.csv(paste(folder_path, "E_coca_124.csv", sep=""), header=TRUE)
ecoca_124 <- ecoca_124[,c(1,2,3,6,18)]

ecoca_209 <- read.csv(paste(folder_path, "E_coca_209.csv", sep=""), header=TRUE)
ecoca_209 <- ecoca_209[,c(1,2,3,6,18)]

ecoca_228 <- read.csv(paste(folder_path, "E_coca_228.csv", sep=""), header=TRUE)
ecoca_228 <- ecoca_228[,c(1,2,3,6,18)]

ehond <- read.csv(paste(folder_path, "E_hond.csv", sep=""), header=TRUE)
ehond <- ehond[,c(1,2,3,6,18)]

#Filter lists to keep only the following: 
#manually annotated pathway genes ("pathway_"), P450 candidates, 2ODD candidates, polyketide candidates

ecoca_48 <- ecoca_48[grep("pathway|20DD|P450|polyketide", ecoca_48$hit, ignore.case = TRUE),]
ecoca_113 <- ecoca_113[grep("pathway|20DD|P450|polyketide", ecoca_113$hit, ignore.case = TRUE),]
ecoca_124 <- ecoca_124[grep("pathway|20DD|P450|polyketide", ecoca_124$hit, ignore.case = TRUE),]
ecoca_209 <- ecoca_209[grep("pathway|20DD|P450|polyketide", ecoca_209$hit, ignore.case = TRUE),]
ecoca_228 <- ecoca_228[grep("pathway|20DD|P450|polyketide", ecoca_228$hit, ignore.case = TRUE),]
ehond <- ehond[grep("pathway|20DD|P450|polyketide", ehond$hit, ignore.case = TRUE),]

#Set e-value of manually curated pathway genes to 0.0.

ecoca_48$evalue[grep("pathway", ecoca_48$hit, ignore.case = TRUE)] <- 0.0
ecoca_113$evalue[grep("pathway", ecoca_113$hit, ignore.case = TRUE)] <- 0.0
ecoca_124$evalue[grep("pathway", ecoca_124$hit, ignore.case = TRUE)] <- 0.0
ecoca_209$evalue[grep("pathway", ecoca_209$hit, ignore.case = TRUE)] <- 0.0
ecoca_228$evalue[grep("pathway", ecoca_228$hit, ignore.case = TRUE)] <- 0.0
ehond$evalue[grep("pathway", ehond$hit, ignore.case = TRUE)] <- 0.0

#Apply evalue cutoff (<=1E-80): remove genes without strong BLAST matches

ecoca_48 <- ecoca_48[ecoca_48$evalue <= 1E-80,]
ecoca_113 <- ecoca_113[ecoca_113$evalue <= 1E-80,]
ecoca_124 <- ecoca_124[ecoca_124$evalue <= 1E-80,]
ecoca_209 <- ecoca_209[ecoca_209$evalue <= 1E-80,]
ecoca_228 <- ecoca_228[ecoca_228$evalue <= 1E-80,]
ehond <- ehond[ehond$evalue <= 1E-50,]

#Sort by evalue, then ensure all hit names are unique (append number for duplicates)

ecoca_48 <- ecoca_48[order(ecoca_48$evalue), ]
ecoca_113 <- ecoca_113[order(ecoca_113$evalue), ]
ecoca_124 <- ecoca_124[order(ecoca_124$evalue), ]
ecoca_209 <- ecoca_209[order(ecoca_209$evalue), ]
ecoca_228 <- ecoca_228[order(ecoca_228$evalue), ]
ehond <- ehond[order(ehond$evalue), ]

ecoca_48$hit <- make.unique(ecoca_48$hit)
ecoca_113$hit <- make.unique(ecoca_113$hit)
ecoca_124$hit <- make.unique(ecoca_124$hit)
ecoca_209$hit <- make.unique(ecoca_209$hit)
ecoca_228$hit <- make.unique(ecoca_228$hit)
ehond$hit <- make.unique(ehond$hit)

#Normalize RPKM values to max RPKM value for each set.

ecoca_48$norm_RPKM <- ecoca_48$RPKM / max(ecoca_48$RPKM)
ecoca_113$norm_RPKM <- ecoca_113$RPKM / max(ecoca_113$RPKM)
ecoca_124$norm_RPKM <- ecoca_124$RPKM / max(ecoca_124$RPKM)
ecoca_209$norm_RPKM <- ecoca_209$RPKM / max(ecoca_209$RPKM)
ecoca_228$norm_RPKM <- ecoca_228$RPKM / max(ecoca_228$RPKM)
ehond$norm_RPKM <- ehond$RPKM / max(ehond$RPKM)

#Sort by norm_RPKM

ecoca_48 <- ecoca_48[order(-ecoca_48$norm_RPKM),]
ecoca_113 <- ecoca_113[order(-ecoca_113$norm_RPKM),]
ecoca_124 <- ecoca_124[order(-ecoca_124$norm_RPKM),]
ecoca_209 <- ecoca_209[order(-ecoca_209$norm_RPKM),]
ecoca_228 <- ecoca_228[order(-ecoca_228$norm_RPKM),]
ehond <- ehond[order(-ehond$norm_RPKM),]

#Compute ratio of norm_RPKM in E. coca vs. E. hondense.
#If gene not found in E. hond, assign it RPKM equal to lowest-expressed gene (to avoid divide-by-0 error).
#In such cases, round ratio to nearest whole number so it is easy to tell which genes are not present in E. hond.

ecoca_48$ratio <- apply(ecoca_48, 1, function(x) 
  {
    if (is.na(match(x[2],ehond$hit))) {round(as.numeric(x[6]) / as.numeric(min(ehond$norm_RPKM)))}
    else {as.numeric(x[6]) / as.numeric(ehond$norm_RPKM[match(x[2],ehond$hit)])}
  })
  
ecoca_113$ratio <- apply(ecoca_113, 1, function(x) 
{
  if (is.na(match(x[2],ehond$hit))) {round(as.numeric(x[6]) / as.numeric(min(ehond$norm_RPKM)))}
  else {as.numeric(x[6]) / as.numeric(ehond$norm_RPKM[match(x[2],ehond$hit)])}
})

ecoca_124$ratio <- apply(ecoca_124, 1, function(x) 
{
  if (is.na(match(x[2],ehond$hit))) {round(as.numeric(x[6]) / as.numeric(min(ehond$norm_RPKM)))}
  else {as.numeric(x[6]) / as.numeric(ehond$norm_RPKM[match(x[2],ehond$hit)])}
})

ecoca_209$ratio <- apply(ecoca_209, 1, function(x) 
{
  if (is.na(match(x[2],ehond$hit))) {round(as.numeric(x[6]) / as.numeric(min(ehond$norm_RPKM)))}
  else {as.numeric(x[6]) / as.numeric(ehond$norm_RPKM[match(x[2],ehond$hit)])}
})

ecoca_228$ratio <- apply(ecoca_228, 1, function(x) 
{
  if (is.na(match(x[2],ehond$hit))) {round(as.numeric(x[6]) / as.numeric(min(ehond$norm_RPKM)))}
  else {as.numeric(x[6]) / as.numeric(ehond$norm_RPKM[match(x[2],ehond$hit)])}
})

#Sort candidates by E. coca / E. hond ratio

ecoca_48 <- ecoca_48[order(-ecoca_48$ratio),]
ecoca_113 <- ecoca_113[order(-ecoca_113$ratio),]
ecoca_124 <- ecoca_124[order(-ecoca_124$ratio),]
ecoca_209 <- ecoca_209[order(-ecoca_209$ratio),]
ecoca_228 <- ecoca_228[order(-ecoca_228$ratio),]

##################################################################
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

#Generate hierarchical clustering / heatmap (use log2 of raw RPKM values)

profile_data <- data.matrix(log2(profile[,-1]))
rownames(profile_data) <- profile$hits

pheatmap(profile_data, cluster_rows=TRUE, cluster_cols=FALSE, scale='row',
         fontsize_row = 3)

#Construct matrix of norm_RPKM values across E. coca and E. hond for each hit 
#(assign norm_RPKM = 0 if gene not found)

profile <- data.frame(hits)

profile$ecoca_48 <- ecoca_48$norm_RPKM[match(hits,ecoca_48$hit)]
profile$ecoca_113 <- ecoca_113$norm_RPKM[match(hits,ecoca_113$hit)]
profile$ecoca_124 <- ecoca_124$norm_RPKM[match(hits,ecoca_124$hit)]
profile$ecoca_209 <- ecoca_209$norm_RPKM[match(hits,ecoca_209$hit)]
profile$ecoca_228 <- ecoca_228$norm_RPKM[match(hits,ecoca_228$hit)]
profile$ehond <- ehond$norm_RPKM[match(hits,ehond$hit)]

profile[is.na(profile)] <- 0

#Generate hierarchical clustering / heatmap

profile_data <- data.matrix(profile[,-1])
rownames(profile_data) <- profile$hits

pheatmap(profile_data, cluster_rows=TRUE, cluster_cols=FALSE, scale='row',
         fontsize_row = 3)
