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
snp10_rec = (snp10 == 2)
glm10_rec = glm(status ~ snp10_rec,family="binomial")
tidy(glm10_rec)
#----------------------------------

#------- question 4 ----------------
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
results <- rep(NA, dim(snpdata)[2]) 
for (i in 1:ncol(snpdata)) {
  snp <- as.numeric(snpdata[,i])
  snp[snp==0] <- NA
  glm_all <- glm(status ~ snp, family = "binomial")
  results[i] <- tidy(glm_all)$statistic[2]
}
mean(results)
min(results)
max(results)
#-----------------------------------

#---- question 5 ----------------------
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]

snpdata = sub.10@.Data
status = subject.support$cc
results <- rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)) {
  snp <- as.numeric(snpdata[,i])
  snp[snp==0] <- NA
  glm <- glm(status ~ snp, family = "binomial")
  results[i] <- glm$coefficients[1] * glm$coefficients[1]
}
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)
slotNames(glm_all)
results_chi <- chi.squared(glm_all)
cor(results_chi, results)
#----------------------------------------

#--- question 6 -------------------------------
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

tstats_obj_pop = rowttests(edata,pdata$population)
fstats_obj_pop = rowFtests(edata,as.factor(pdata$population))

tstats_obj_study = rowttests(edata,pdata$study)
fstats_obj_study = rowFtests(edata,as.factor(pdata$study))

cor(tstats_obj_pop$statistic, fstats_obj_pop$statistic)
#[1] -0.9502577
cor(tstats_obj_pop$p.value, fstats_obj_pop$p.value)
#[1] 1
cor(tstats_obj_study$statistic, fstats_obj_study$statistic)
#[1] -0.9502577
cor(tstats_obj_study$p.value, fstats_obj_study$p.value)
#[1] 1
#----------------------------------------------

#---- question 7 -------------------------------------------
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)
library(MASS)
library(DESeq2)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)

edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ pdata$study)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
cor(limma_pvals, result_nb$pvalue)
#-----------------------------------------------------------
