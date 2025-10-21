## Proportions of ambiguous sites
library(data.table)
m1A=read.csv("m6Anet_m1A/data.indiv_proba.csv")
m6A=read.csv("m6Anet/data.indiv_proba.csv")
AI=read.csv("m6Anet_edi/data.indiv_proba.csv")
colnames(m1A)[4]="m1A_ratio"
colnames(m6A)[4]="m6A_ratio"
colnames(AI)[4]="AI_ratio"
merged_data <- merge(m1A, m6A, by = c("transcript_id", "transcript_position", "read_index"))
merged_data <- merge(merged_data, AI, by = c("transcript_id", "transcript_position", "read_index"))
m6Anet <- merged_data[, c("transcript_id", "transcript_position", "read_index", "m1A_ratio", "m6A_ratio", "AI_ratio")]

m1A=read.table("TandemMod_m1A/predict_m1A.tsv",header=F,sep="\t")
AI=read.table("TandemMod_edi/predict_AtoI.tsv",header=F,sep="\t")
m6A=read.table("TandemMod/predict_m6A.tsv",header=F,sep="\t")
colnames(m1A)[6]="m1A_ratio"
colnames(m6A)[6]="m6A_ratio"
colnames(AI)[6]="AI_ratio"
merged_data <- merge(m1A, m6A, by = c("V1", "V2", "V4"))
merged_data <- merge(merged_data, AI, by = c("V1", "V2", "V4"))
Tandem <- merged_data[, c("m1A_ratio", "m6A_ratio", "AI_ratio")]

m1A=read.table("SingleMod_m1A/prediction_m1A.txt",header=F,sep="\t")
AI=read.table("SingleMod_edi/prediction_AI.txt",header=F,sep="\t")
m6A=read.table("SingleMod_m6A/prediction_m6A.txt",header=F,sep="\t")
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

##signal-derived features between true positive and false positive sites identified by TandemMod
merged_data <- read.delim("TandemMod_featrue.tsv", header = TRUE)
AI_high <- merged_data
m6A_high <- merged_data
m1A_high <- merged_data 
AI_site=read.table("AtoI_groundtruth.bed",sep="\t",header=F)
colnames(AI_site)=c("read_id","transcript_id","transcript_position")
AI_high$group <- ifelse(
  paste(AI_high$transcript_id, AI_high$transcript_position, sep = "_") %in%
    paste(AI_site$transcript_id, AI_site$transcript_position, sep = "_"),
  "groundtruth",
  "non-groundtruth"
)
m1A_site=read.table("m1A_groundtruth.bed",sep="\t",header=F)
m6A_site=read.table("m6A_groundtruth.bed",sep="\t",header=F)
colnames(m1A_site)=c("read_id","transcript_id","transcript_position")
m1A_high$group <- ifelse(
  paste(m1A_high$transcript_id, m1A_high$transcript_position, sep = "_") %in%
    paste(m1A_site$transcript_id, m1A_site$transcript_position, sep = "_"),
  "groundtruth",
  "non-groundtruth"
)
colnames(m6A_site)=c("read_id","transcript_id","transcript_position")
m6A_high$group <- ifelse(
  paste(m6A_high$transcript_id, m6A_high$transcript_position, sep = "_") %in%
    paste(m6A_site$transcript_id, m6A_site$transcript_position, sep = "_"),
  "groundtruth",
  "non-groundtruth"
)
library(tidyverse)
library(ggplot2)
combined_data <- bind_rows(
  m6A_high %>% mutate(type = "m6A"),
  m1A_high %>% mutate(type = "m1A"),
  AI_high %>% mutate(type = "AtoI")
)
df_long <- combined_data %>%
  select(type, group, std_1, std_2, std_3, std_4, std_5) %>%
  pivot_longer(cols = starts_with("std_"),
               names_to = "std_index",
               values_to = "value") %>%
  mutate(std_num = as.numeric(str_extract(std_index, "\\d+")))
df_long <- combined_data %>%
  select(type, group, mean_1, mean_2, mean_3, mean_4, mean_5) %>%
  pivot_longer(cols = starts_with("mean_"),
               names_to = "mean_index",
               values_to = "value") %>%
  mutate(std_num = as.numeric(str_extract(mean_index, "\\d+")))

df_long <- combined_data %>%
  select(type, group, length_1, length_2, length_3, length_4, length_5) %>%
  pivot_longer(cols = starts_with("length_"),
               names_to = "length_index",
               values_to = "value") %>%
  mutate(std_num = as.numeric(str_extract(length_index, "\\d+")))

df_long <- combined_data %>%
  select(type, group, kmer_base_quality_1, kmer_base_quality_2, kmer_base_quality_3, kmer_base_quality_4, kmer_base_quality_5) %>%
  pivot_longer(cols = starts_with("kmer_base_quality_"),
               names_to = "kmer_base_quality_index",
               values_to = "value") %>%
  mutate(std_num = as.numeric(str_extract(kmer_base_quality_index, "\\d+")))

summary_df <- df_long %>%
  group_by(type, group, std_num) %>%
  summarise(
    q1 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(summary_df, aes(x = std_num, color = group)) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = 0,
    position = position_nudge(x = ifelse(summary_df$group == "groundtruth", -0.1, 0.1)),
    size = 0.3
  ) +
  geom_errorbar(
    aes(ymin = q1, ymax = q1),
    width = 0.2,
    position = position_nudge(x = ifelse(summary_df$group == "groundtruth", -0.1, 0.1)),
    size = 0.3
  ) +
  geom_errorbar(
    aes(ymin = q3, ymax = q3),
    width = 0.2,
    position = position_nudge(x = ifelse(summary_df$group == "groundtruth", -0.1, 0.1)),
    size = 0.3
  ) +
  geom_point(
    aes(y = median),
    position = position_nudge(x = ifelse(summary_df$group == "groundtruth", -0.1, 0.1)),
    size = 1.5
  ) +
  scale_x_continuous(
    breaks = 1:5,
    labels = 1:5,
    limits = c(0.5, 5.5),
    expand = expansion(mult = 0)
  ) +
  labs(x = "k-mer base", y = "std") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid(. ~ type, scales = "free_x", space = "free_x")  
ggsave("TandemMod_std.pdf",p,width=4,height=2.5)