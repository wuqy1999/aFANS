library(FeaAug)
#the package "FeaAug" can be downloaded from https://github.com/ShanLu92/FeaAug
library(MASS)
library(ROCR)
library(cutpointr)
library(tidyverse)
library(DIRECT)
source("aFANS_Functions.R",encoding = "utf-8")

# load the sample dataset
load("example_data.RData")

#apply PL50 test to the dataset
PL50_var_judge(data11,yname="species",y=c(1,0),alpha = 0.05,iter = 1000)
##$same[1]  1  2  3  4  5  6  7  8  9 10
##$diff[1] 11 12 13 14 15 16 17 18 19 20

#test the efficiency of LR
lr <- lr_rate(
  dataset = data11,
  iter = 10,
  subtrain = 300,
  subtest = 300
)
#test the efficiency of FANS
fans <-
  fans_lr_rate(
    dataset = data11,
    iter = 10,
    subtrain = 300,
    subtest = 300,
    L = 5
  )
#test the efficiency of aFANS
afans <-
  afans_lr_rate(
    dataset = data11,
    iter = 10,
    subtrain = 300,
    subtest = 300,
    L = 5,
    alpha_var = 0.05,
    var_iter = 1000
  )
#test the efficiency of FANS2
fans2 <-
  fans2_lr_rate(
    dataset = data11,
    iter = 10,
    subtrain = 300,
    subtest = 300,
    L = 5
  )
#print the accuracy, false positive rate, false negative rate, and AUC of each iteration(aFANS).
(afans$rate)
#print the time of each iteration
(afans$time)
