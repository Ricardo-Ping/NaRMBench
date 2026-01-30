pacman::p_load(data.table,tidyverse,caret,stringr,keras,tensorflow, optparse, multiROC,pROC,PRROC)

rm(list=ls())
gc()

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file name", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of cores allocated", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be supplied (-i).", call.=FALSE)
}
getdinodir <- function(){
    commandArgs() %>%
       tibble::enframe(name=NULL) %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value) %>%
       word(., start=1, end=-3, sep="/")
}
dinodir <- getdinodir()

Rcpp::sourceCpp(paste0(dinodir,"/code/all_functions.cpp"))

chrmapping <- readRDS("chrmapping.rds")

model1=load_model_hdf5(paste0(dinodir, "/model/best_pos5_mix_3class_resnet_1992.h5"))
model2=load_model_hdf5(paste0(dinodir, "/model/best_pos5_mix_3c_1vs1_3010_resnet_10.h5"))
model3=load_model_hdf5(paste0(dinodir, "/model/best_regression_mixnHEK_morefts_16384_1024_b1024_init_650k_XHe0_Mus_asin06.h5"))

load(opt$input)

agggrp <- word(opt$input,sep = "\\.", 1)
  
k_set_learning_phase(0)

pred <- predict(model1,test$x[,,c(1:43),])
pred <- as.data.frame(pred)

pred_test <- colnames(pred[,c(1:3)])[apply(pred[,c(1:3)],1,which.max)]
pred_test <- gsub("V1",0,pred_test)
pred_test <- gsub("V2",1,pred_test)
pred_test <- gsub("V3",2,pred_test)
pred_test <- as.numeric(pred_test)

c01=which(pred_test!=2)
pred$n2=1-pred$V3
pred$c1=0
pred$c2=0
  
pred2 <- predict(model2,test$x[c01,,c(1:43),])
pred2 <- as.data.frame(pred2)
  
pred[c01,c(5,6)] <- pred2
pred[c01,"V1"] <- pred[c01,"c1"]*pred[c01,"n2"]
pred[c01,"V2"] <- pred[c01,"c2"]*pred[c01,"n2"]

pred.2model <- data.frame(pred[,c(1:3)], 'ref'=test$y, 'id'=test$info[,"id"], 'cov'=test$info[,"cov"]) %>% 
  group_by(id, cov, ref) %>% 
  dplyr::summarise('0'=mean(V1),'1'=mean(V2),'2'=mean(V3))

pred.2model$pred <- colnames(pred.2model[,c(4:6)])[apply(pred.2model[,c(4:6)],1,which.max)]
pred.2model <- pred.2model %>% 
  separate(id,into=c("chr_str","position"), sep=":", remove = F) %>% 
  mutate(strand = str_sub(chr_str,-1),
         strand = ifelse(strand == 1, "+", "-"),
         chr = str_sub(chr_str, 1, -2),
         cov = as.numeric(cov))

pred.2model$contig=replace_list(chrmapping$chromid,chrmapping$chroms,pred.2model$chr)


getROC <- function(df) {
  # Select relevant columns for "S0" and "S1" classes
  #y_ref <- to_categorical(df$ref, 2) # Only two classes: S0 and S1
  perform.df <- data.frame(df$ref, df[, 6])
  
  colnames(perform.df) <- c("S1_true","S1_pred_n")
  roc_curve <- roc(response = perform.df$S1_true, predictor = perform.df$S1_pred_n)
  roc_df <- data.frame(
    FPR = 1 - roc_curve$specificities,  # False Positive Rate
    TPR = roc_curve$sensitivities       # True Positive Rate
  )
  auc_value <- auc(roc_curve)
  roc_df$AUC <- auc_value
  
  # Append (0, 0) to the beginning if not already present
  if (!(0 %in% roc_df$FPR & 0 %in% roc_df$TPR)) {
    roc_df <- rbind(data.frame(FPR = 0, TPR = 0, AUC = auc_value), roc_df)
  }

  # Append (1, 1) to the end if not already present
  if (!(1 %in% roc_df$FPR & 1 %in% roc_df$TPR)) {
    roc_df <- rbind(roc_df, data.frame(FPR = 1, TPR = 1, AUC = auc_value))
  }

  return(roc_df)
}

# Process one-model and two-model data
#roc.1model <- getROC(pred.1model) %>% 
  #mutate(model = paste0("One model (AUC:", round(AUC[1], 3), ")"))
#roc.2model <- getROC(pred.2model) %>% 
  #mutate(model = paste0("models (AUC:", round(AUC[1], 3), ")"))
#roc.all <- rbind(roc.1model, roc.2model)

roc.2model <- getROC(pred.2model)
auc_value <- unique(roc.2model$AUC)  # Ensure AUC is unique
if (length(auc_value) != 1) {
  stop("AUC values in the ROC data frame are not consistent.")
}

roc.2model <- roc.2model %>% 
  mutate(model = paste0("models (AUC:", round(auc_value, 3), ")"))


# Create the ROC plot
roc = ggplot(data = roc.2model, aes(x = FPR, y = TPR, colour = model)) +
  geom_line(size = 1) +
  scale_color_manual(values=c("#A681BA")) + 
  theme_bw() +
  ylim(0, 1) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = "dashed") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(legend.position = c(0.7, 0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(2, "line"),
        axis.text.y = element_text(size = 20, color = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 20, color = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title = element_text(size = 25, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.ticks.length = unit(1, "line"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save the plot
ggsave("ROC_2class_plot.pdf", plot = roc, width = 8, height = 6)

getPR <- function(df) {

  perform.df <- data.frame(df$ref, df[, 6])
  colnames(perform.df) <- c("true", "pred")

  pr_curve <- pr.curve(scores.class0 = perform.df$pred[perform.df$true == 1],
                       weights.class0 = perform.df$true,
                       curve = TRUE)

  pr_df <- data.frame(
    recall = pr_curve$curve[, 1],
    precision = pr_curve$curve[, 2]
  )
  pr_df$AUC <- pr_curve$auc.integral

  return(pr_df)
}

# Calculate PR curve for the given model
pr.2model <- getPR(pred.2model)
auc_value <- unique(pr.2model$AUC)  # Ensure AUC is unique
if (length(auc_value) != 1) {
  stop("Multiple inconsistent AUC values found in the PR data frame.")
}

pr.2model <- pr.2model %>% 
  mutate(model = paste0("models (AUC:", round(auc_value, 3), ")"))

# Create and save the PR curve plot
pr <- ggplot(data = pr.2model, aes(x = recall, y = precision, colour = model)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#A681BA")) + 
  theme_bw() +
  ylim(0, 1) +
  xlab("Recall") +
  ylab("Precision") +
  theme(legend.position = c(0.3, 0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width = unit(2, "line"),
        legend.key.height = unit(2, "line"),
        axis.text.y = element_text(size = 20, color = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 20, color = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title = element_text(size = 25, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.ticks.length = unit(1, "line"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("PR_2class_plot.pdf", plot = pr, width = 8, height = 6)


#####Regression model to predict editing rate####
if (FALSE) {
pred.2modelcov40 <- pred.2model %>% filter(cov >= 40)
id=pred.2modelcov40$id[pred.2modelcov40$pred %in% c("1")]
ind=which(test$info[,"id"] %in% id)
pred3 <- predict(model3,0.01*test$x[ind,,,])

pred.regression <- data.frame('id'=(test$info[ind,"id"]),
                              'cov'=test$info[ind,"cov"],
                              'ref'=test$y[ind],
                              'rate'=test$y2[ind],
                              'pred.rate'=pred3[,1]) %>% 
  group_by(id,cov,ref, rate) %>% 
  dplyr::summarise(pred.rate=mean(pred.rate)) %>% 
  ungroup()

#Combine classification and regression model together
pred.all <- left_join(pred.2model, pred.regression[,c(1,4,5)],by="id") %>% 
  filter(cov >= 20) %>% 
  mutate(pred.rate = sin(pred.rate**(5/3)))
}

#class0 <- pred.all[pred.all$pred == 0,c("contig","position","strand","cov","0","1","2","ref","pred")]
#class1 <- pred.all[pred.all$pred == 1,c("contig","position","strand","cov","0","1","2","ref","pred","rate","pred.rate")]
#class2 <- pred.all[pred.all$pred == 2,c("contig","position","strand","cov","0","1","2","ref","pred")]

#fwrite(class0, paste0(agggrp, ".output_prediction_CNN_class0.txt"),sep = "\t")
#fwrite(class1, paste0(agggrp, ".output_prediction_CNN_class1.txt"),sep = "\t")
fwrite(pred.2model, paste0(agggrp, ".output_prediction_CNN_class2.txt"),sep = "\t")
