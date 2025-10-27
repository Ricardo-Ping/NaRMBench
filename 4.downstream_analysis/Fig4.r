library(dplyr)
library(ggplot2)
library(tidyr)
library(metaplot)
library(ggprism)
library(ggplot2)

##Distribution of predicted modification sites across transcript regions
bed_files <- list.files(path = ".", pattern = "*.txt$", full.names = TRUE)
#GroundTruth: ../ground_truth_sites/
Group=c("Dorado-RNA004","Epinano-RNA004","GroundTruth","m6Anet-RNA004","SingleMod-RNA004","TandemMod-RNA004","Xron-RNA004")
colors = c("#6affb9", "#9fcc62","black", "#f6c365", "#2e3792", "#8e5aa2", "#ef1fff")
feature_anno=prepare_features(gtf_file="Homo_sapiens.GRCh38.112.gtf.gz")
meta=metaPlot(bed_file=bed_files,group_names=Group,features_anno=feature_anno,scale_region=T)
distru=meta$data
p <- ggplot(distru, aes(x = rel_pos, color = group_names)) +
  geom_density(size = 0.7) +  
  geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "grey40", size = 1.2) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  
  scale_x_continuous(breaks = c(0.5, 1.5, 2.5),  
                     labels = c("5'UTR", "CDS", "3'UTR"),
                     limits = c(0, 3)) +  
  scale_color_manual(values = setNames(colors, Group)) +  
  theme(axis.line = element_line(arrow = arrow(length = unit(0.25, "cm"), type = "closed")),  
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(face = "bold")) +
  ylab("Density") +  
  xlab("") +
  theme_classic(base_size = 12)
ggsave("m6A_distribution.pdf",p,height=7,width=10)

##WT vs. KO (or unmodified IVT) comparison metrics

# Define tool thresholds
tool_thresholds <- list(
    'Xron-RNA004'= 0.01,
    'Dorado-RNA004'= 0.2,
    'SingleMod-RNA004'= 0.17,
    'TandemMod-RNA004'=0.68,
    'm6Anet-RNA004'=0.75,
    'EpiNano-RNA004'=0.48
)

tools <- c("Xron-RNA004", "Dorado-RNA004","SingleMod-RNA004","TandemMod-RNA004","m6Anet-RNA004",'EpiNano-RNA004')

# Initialize data storage
output_data <- data.frame(
  Tool = character(),
  WT_KO_WT = numeric(),
  WT_gt_KO_WT_KO = numeric(),
  WT_count = numeric(),
  KO_count = numeric(),
  stringsAsFactors = FALSE
)

# Process each tool
for (tool in tools) {
  ko_file <- paste0("./ko/",tool, "_m6A_with_labels_converted.txt")
  wt_file <- paste0("./", tool, "_m6A_with_labels_converted.txt")
  data <- read.table(wt_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  wt_data <- data %>%
  separate(V1, into = c("chr", "start"), sep = "_") %>%
  mutate( 
    end = start  
  ) %>%
  select(chr, start, end, status = V2, value = V3)
  
  data <- read.table(ko_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  ko_data <- data %>%
  separate(V1, into = c("chr", "start"), sep = "_") %>%
  mutate(
    end = start  
  ) %>%
  select(chr, start, end, status = V2, value = V3)
  
  colnames(ko_data) <- c('chr', 'start', 'end', 'status', 'value')
  colnames(wt_data) <- c('chr', 'start', 'end', 'status', 'value')
  threshold <- tool_thresholds[[tool]]
  
  wt_filtered <- wt_data %>% filter(value >= threshold)
  ko_filtered <- ko_data %>% filter(value >= threshold)
  
  wt_ko_common <- inner_join(wt_filtered, ko_filtered, 
                            by = c('chr', 'start', 'end'),
                            suffix = c('_wt', '_ko'))
  
  wt_count <- nrow(wt_filtered)
  wt_ko_count <- nrow(wt_ko_common)
  ko_count <- nrow(ko_filtered)
  
  wt_ko_ratio <- ifelse(wt_count > 0, wt_ko_count / wt_count, 0)
  
  if (tool == "ELIGOS") {
    wt_gt_ko_ratio <- ifelse(wt_ko_count > 0, 
                            sum(wt_ko_common$value_wt < wt_ko_common$value_ko) / wt_ko_count, 
                            0)
  } else {
    wt_gt_ko_ratio <- ifelse(wt_ko_count > 0, 
                            sum(wt_ko_common$value_wt > wt_ko_common$value_ko) / wt_ko_count, 
                            0)
  }
  
  output_data <- rbind(output_data, data.frame(
    Tool = tool,
    WT_KO_WT = wt_ko_ratio,
    WT_gt_KO_WT_KO = wt_gt_ko_ratio,
    WT_count = wt_count,
    KO_count = ko_count
  ))
}

# Reorder tools by WT_KO_WT ratio descending
output_data <- output_data %>% arrange(desc(WT_KO_WT))
tools_ordered <- output_data$Tool

# Convert to long format for plotting
plot_data <- output_data %>%
  select(Tool, WT_KO_WT, WT_gt_KO_WT_KO) %>%
  pivot_longer(cols = -Tool, names_to = "Type", values_to = "Value")
#if(1==2){
# Plot the ratios
ggplot(plot_data, aes(x = factor(Tool, levels = tools_ordered), y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("WT_KO_WT" = "lightblue", "WT_gt_KO_WT_KO" = "orange"),
                    labels = c("(WT>KO)/WT&KO","(WT&KO)/WT")) +
  labs(y = "Ratio",x=NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Comparison of WT&KO/WT and WT>KO/WT&KO Ratios")

ggsave("m6A_difference_barplot.pdf", width = 8, height = 5)

ratio_data <- output_data %>%
  mutate(ratio = KO_count / WT_count,
         ratio_log2 = log2(ratio)) 

ggplot(ratio_data, aes(x = factor(Tool, levels = tools_ordered), y = ratio_log2, fill = ratio_log2 > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +  
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(y = "log2(KO / WT)", x = NULL, fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("m6A_count_barplot.pdf", width = 8, height = 5)
