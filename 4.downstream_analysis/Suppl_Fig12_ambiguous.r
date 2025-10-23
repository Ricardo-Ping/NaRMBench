## Proportions of ambiguous sites
library(data.table)
m1A=read.csv("./example_data/m6Anet_m1A_read_proba.csv")
m6A=read.csv("./example_data/m6Anet_m6A_read_proba.csv")
AI=read.csv("./example_data/m6Anet_edi_read_proba.csv")
colnames(m1A)[4]="m1A_ratio"
colnames(m6A)[4]="m6A_ratio"
colnames(AI)[4]="AI_ratio"
merged_data <- merge(m1A, m6A, by = c("transcript_id", "transcript_position", "read_index"))
merged_data <- merge(merged_data, AI, by = c("transcript_id", "transcript_position", "read_index"))
m6Anet <- merged_data[, c("transcript_id", "transcript_position", "read_index", "m1A_ratio", "m6A_ratio", "AI_ratio")]

m1A=read.table("./example_data/TandemMod_m1A_read_proba.tsv",header=F,sep="\t")
AI=read.table("./example_data/TandemMod_edi_read_proba.tsv",header=F,sep="\t")
m6A=read.table("./example_data/TandemMod_m6A_read_proba.tsv",header=F,sep="\t")
colnames(m1A)[6]="m1A_ratio"
colnames(m6A)[6]="m6A_ratio"
colnames(AI)[6]="AI_ratio"
merged_data <- merge(m1A, m6A, by = c("V1", "V2", "V4"))
merged_data <- merge(merged_data, AI, by = c("V1", "V2", "V4"))
Tandem <- merged_data[, c("m1A_ratio", "m6A_ratio", "AI_ratio")]

m1A=read.table("./example_data/SingleMod_m1A_read_proba.txt",header=F,sep="\t")
AI=read.table("./example_data/SingleMod_edi_read_proba.txt",header=F,sep="\t")
m6A=read.table("./example_data/SingleMod_m6A_read_proba.txt",header=F,sep="\t")
colnames(m1A)[2]="m1A_ratio"
colnames(m6A)[2]="m6A_ratio"
colnames(AI)[2]="AI_ratio"
merged_data <- merge(m1A, m6A, by = c("V1"))
merged_data <- merge(merged_data, AI, by = c("V1"))
singlemod <- merged_data[, c("m1A_ratio", "m6A_ratio", "AI_ratio")]

library(UpSetR)
sapply(Tandem[, c("m1A_ratio", "m6A_ratio", "AI_ratio")], function(x) quantile(x, 0.9))
data <- Tandem
data$m1A <- data$m1A_ratio > 0.9
data$m6A <- data$m6A_ratio > 0.9
data$AtoI <- data$AI_ratio > 0.9
data <- data.frame(
  m1A = as.numeric(data$m1A),
  m6A = as.numeric(data$m6A),
  AtoI = as.numeric(data$AtoI)
)
pdf("tandem.pdf", height = 5, width = 8)
upset(data, sets = c("m6A","m1A","AtoI"), keep.order = TRUE,order.by = c("degree"),point.size = 1.5,line.size = 0.5,sets.bar.color = "grey",main.bar.color = 'grey',text.scale = c(1.3, 1.3, 1, 1, 2, 1.8),mainbar.y.label = "Modificaton Intersections", sets.x.label = "Size per Mod")
dev.off()
data <- m6Anet
data$m1A <- data$m1A_ratio > 0.9
data$m6A <- data$m6A_ratio > 0.9
data$AtoI <- data$AI_ratio > 0.9
data <- data.frame(
  m1A = as.numeric(data$m1A),
  m6A = as.numeric(data$m6A),
  AtoI = as.numeric(data$AtoI)
)
pdf("m6Anet.pdf", height = 5, width = 8)
upset(data, sets = c("m6A","m1A","AtoI"), keep.order = TRUE,order.by = c("degree"),point.size = 1.5,line.size = 0.5,sets.bar.color = "grey",main.bar.color = 'grey',text.scale = c(1.3, 1.3, 1, 1, 2, 1.8),mainbar.y.label = "Modificaton Intersections", sets.x.label = "Size per Mod")
dev.off()
data <- singlemod
data$m1A <- data$m1A_ratio > 0.9
data$m6A <- data$m6A_ratio > 0.9
data$AtoI <- data$AI_ratio > 0.9
data <- data.frame(
  m1A = as.numeric(data$m1A),
  m6A = as.numeric(data$m6A),
  AtoI = as.numeric(data$AtoI)
)
pdf("singlemod.pdf", height = 5, width = 8)
upset(data, sets = c("m6A","m1A","AtoI"), keep.order = TRUE,order.by = c("degree"),point.size = 1.5,line.size = 0.5,sets.bar.color = "grey",main.bar.color = 'grey',text.scale = c(1.3, 1.3, 1, 1, 2, 1.8),mainbar.y.label = "Modificaton Intersections", sets.x.label = "Size per Mod")
dev.off()

