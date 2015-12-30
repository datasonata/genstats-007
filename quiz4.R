library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)
library(limma)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]


edata_sub <- edata[,c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19)]
pdata_bm_nona <- pdata_bm[!is.na(pdata_bm$age)==TRUE,]

mod_adj = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata,mod_adj)
limma_stats <- fit_limma$coefficients
nrow(limma_stats[as.integer(limma_stats[,1]) >= 0.05,])

de = DESeqDataSetFromMatrix(edata, pdata_bot, ~strain)
de_fit = DESeq(de)
de_results = results(de_fit)
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes)
names(genes) = rownames(edata)
genes = genes[not_na]

pwf=nullp(genes,"mm9","ensGene")
head(pwf)
GO.wall=goseq(pwf,"mm9","ensGene")
head(GO.wall[ order(-GO.wall$over_represented_pvalue), ])
