library(writexl)

#Load in analyzed dataset .RData file

load("~/Desktop/Prashanth/Labwork/Erythroxylaceae project/Sequences/Batch 2 - February 2021/Complete_dataset/Ecoca_diffexp_2021-05-02.RData")

#Load in list of candidates from heatmaps

candidates <- read.table("~/Desktop/Prashanth/Labwork/Erythroxylaceae project/Sequences/Batch 2 - February 2021/Complete_dataset/Ecoca_diffexp_heatmap_candidates_v1.txt", sep = '\n')

#Make and fill dataframe to hold hit names and corresponding locus IDs for each dataset

locusIDs <- data.frame(hit_name = candidates$V1)

locusIDs$ecoca_48 <- ecoca_48$ID[match(locusIDs$hit_name, ecoca_48$hit)]
locusIDs$ecoca_113 <- ecoca_113$ID[match(locusIDs$hit_name, ecoca_113$hit)]
locusIDs$ecoca_124 <- ecoca_124$ID[match(locusIDs$hit_name, ecoca_124$hit)]
locusIDs$ecoca_209 <- ecoca_209$ID[match(locusIDs$hit_name, ecoca_209$hit)]
locusIDs$ecoca_228 <- ecoca_228$ID[match(locusIDs$hit_name, ecoca_228$hit)]

#Write table to excel file

write_xlsx(locusIDs, "~/Desktop/Prashanth/Labwork/Erythroxylaceae project/Sequences/Batch 2 - February 2021/Complete_dataset/Ecoca_diffexp_heatmap_candidates_v1.xlsx")
