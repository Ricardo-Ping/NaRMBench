library(fmsb)
library(dplyr)
library(tibble)
#radar data : Supplementary Table 4 
data=read.csv("m6A_redar.csv",check.names = FALSE)
data_normalized <- data
cols_invert_1_minus <- c("Motif bias", "Depth bias", "Modification level bias")
for (col in cols_invert_1_minus) {
    data_normalized[[col]] <- 1 - data_normalized[[col]]
}
cols_invert_max_minus <- c("Speed", "Memory efficiency","Differen between KO and WT")
for (col in cols_invert_max_minus) {
    data_normalized[[col]] <- max(data_normalized[[col]]) - data_normalized[[col]]
}
for (col in cols_invert_max_minus) {
    data_normalized[[col]] <- (data_normalized[[col]] - min(data_normalized[[col]])) /
                              (max(data_normalized[[col]]) - min(data_normalized[[col]]))
}
max_values <- apply(data_normalized[, -1], 2, max)  
min_values <- apply(data_normalized[, -1], 2, min)  
max_row <- c("Max", rep(1, length(max_values)))  
min_row <- c("Min", rep(0, length(min_values)))
data_radar <- rbind(max_row, min_row, data_normalized)
data_radar <- column_to_rownames(data_radar, var = "Tool")
data_radar_numeric <- as.data.frame(apply(data_radar, c(1, 2), function(x) {
  num <- suppressWarnings(as.numeric(x))
  if (is.na(num) && !is.na(x) && x != "NA" && x != "<NA>") {
    return(x) 
  } else {
    return(num)
  }
}))
colnames(data_radar_numeric)=c("AUROC","AUPRC","Distribution similarity with ground-truth","Differen between KO and WT","Conserved motif bias","Sequencing depth bias","Motification level bias","Speed","Memory efficiency")
pdf("m6A_radar.pdf",width=3,height=3)
radarchart(
  data_radar_numeric,
  pcol = tool_colors, plwd = 0.5, plty = 1,pty=32,
  cglcol = "grey90", cglty = 1, cglwd = 0.3,
  axislabcol = "grey",
  vlcex = 0.2, vlabels = colnames(data_radar_numeric),
  caxislabels = c("","","","",""),axistype = 0)

dev.off()
