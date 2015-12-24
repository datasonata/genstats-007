library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(UsingR)
library(broom)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

edata = as.matrix(edata)
lm3 = lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
tidy(lm3)

#- question 2 --------------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)
edata = edata - rowMeans(edata)
set.seed(333)
kmeans2 = kmeans(t(edata), centers=2)
svd2 = svd(edata)
cor(kmeans2$cluster, svd2$v)[1]
#----------------------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod,t(edata))
names(fit)
dim(fit$residuals)
dim(fit$effects)
dim(fit$coefficients)
head(fit$effects)

str(fit$effects)

#-- question 7 -----------------------------------------------------------
library(devtools)
library(Biobase)
library(limma)
library(edge)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
#edata = log2(as.matrix(edata) + 1)
#edata = edata[rowMeans(edata) > 10, ]
edata = edata[,c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19)]
mod_adj = model.matrix(~ pdata_bm$age)
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)
fit_limma$coefficients[1000,]
fit_adj$coefficients[,1]
plot(edata[1000,],breaks=100,col=2,xlab="Intercept")

#- -- Question 10 -------------------------------------------------------------
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
#edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) >= 1,]
edata = log2(edata + 1)
edata = edata[,c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19)]
mod = model.matrix(~age,data=pdata_bm)
mod0 = model.matrix(~1, data=pdata_bm)
set.seed(33353)
sva1 = sva(edata,mod,n.sv=2)
