library(dplyr)
library(ggplot2)
library(reshape2)
library(magrittr)
library(ggpubr)
#library(ggpmisc)
library(reshape2)
library(tidyverse)
library(viridis)
library(lme4)
library(lmerTest)
library(Rmisc)
library(data.table)
#library(coin)
library(rstatix)

baseDir='/data/Wangjing/projects/priority/'
subjs=sprintf("%02d", 2:30)


#day2_D1 test
for (s in subjs){
  sub<-paste('sub-',s,sep='')
  if (s==subjs[1]) {
    d2df<-read.csv(paste(baseDir, 'rawdata/', sub, '/ses-02/func/', sub, '_ses-02_task-testD2_events.tsv', sep = ''), header = T, sep = '\t')
    d2df<-subset(d2df, selec=c(rep, condition, word, picID,pic, choice1, choice2, choice3, choice4,pic_choice, choice_RT, hits, cat_group, exp_phase))
    d2df$sub<-sub
  }
  else{
    temp<-read.csv(paste(baseDir, 'rawdata/', sub, '/ses-02/func/', sub, '_ses-02_task-testD2_events.tsv', sep = ''), header = T, sep = '\t')
    temp<-subset(temp, selec=c(rep, condition, word, picID,pic, choice1, choice2, choice3, choice4,pic_choice, choice_RT, hits, cat_group, exp_phase))
    temp$sub<-sub
    d2df<-rbind(d2df, temp)
  }
}  



####d2####
d2df$tested<-ifelse(str_detect(d2df$condition, '_T'), 'test','restudy')

d2_acc_sub_rep<-reshape2::dcast(d2df, sub~rep, value.var = "hits", mean)
d2_RT_sub_rep<-reshape2::dcast(subset(d2df, d2df$hits==1), sub~rep, value.var = "choice_RT", median)

d2_acc_sub_allcond<-reshape2::dcast(d2df, sub~rep+tested, value.var = "hits", mean)
d2_RT_sub_allcond<-reshape2::dcast(subset(d2df, d2df$hits==1), sub~rep+tested, value.var = "choice_RT", median)

view(d2_acc_sub_rep)
mean(d2_acc_sub_rep$"1")
