pacman::p_load(data.table,tidyverse,caret,stringr,keras,tensorflow, optparse, multiROC)

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

class2 <- pred.2model[,c("contig","position","strand","cov","1","pred")]

fwrite(class2, paste0(agggrp, ".prediction_all.txt"),sep = "\t")
