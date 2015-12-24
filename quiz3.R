library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)

#---- question 1 -----------------
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA

lm3 = lm(status ~ snp3)
tidy(lm3)
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)
#----------------------------------

#----- question 3 -----------------
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
snp10 = as.numeric(snpdata[,10])
snp1_dom = (snp1 == 1)
glm1_dom = glm(status ~ snp1_dom,family="binomial")
tidy(glm1_dom)
#----------------------------------

#------- question 4 ----------------
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
snpdata[snpdata==0] = NA
status = subject.support$cc
glm = glm(status ~ as.matrix(snpdata))
slotNames(glm)
glm_all
#-----------------------------------
