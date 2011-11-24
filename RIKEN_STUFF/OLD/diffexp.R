#!/usr/bin/R

library(limma)

rep = read.delim(file='overlap_rna_rep_PYR.txt')
norm = read.delim(file='norm.txt')
coord = read.delim(file='coord.txt')
target = read.delim(file='target.txt')

target$chem_time = paste(target$chem,target$time,sep='_')

exp.rep = norm[norm$id %in% rep$id,]

design = model.matrix(~0+target$chem_time)
colnames(design) = gsub('target$chem_time','',colnames(design),fixed=T)
eset = exp.rep[,-1]
rownames(eset) = exp.rep$id
fit = lmFit(eset, design)

# TSA_2d TSA_2h VEHtsa_2d VEHtsa_2h VEHvpa_2d VEHvpa_2h VPA_2d VPA_2h

# CHANGE BETWEEN TIME
contrast.matrix = makeContrasts(TSA_2h-VEHtsa_2h,levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit3 = eBayes(fit2)
topTable(fit3,lfc=log2(1.5),p.value=0.05,adjust="BH",number=Inf)

contrast.matrix = makeContrasts(TSA_2d-VEHtsa_2d,levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit3 = eBayes(fit2)
topTable(fit3,lfc=log2(1.5),p.value=0.05,adjust="BH",number=Inf)

contrast.matrix = makeContrasts(VPA_2h-VEHvpa_2h,levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit3 = eBayes(fit2)
topTable(fit3,lfc=log2(1.5),p.value=0.05,adjust="BH",number=Inf)

contrast.matrix = makeContrasts(VPA_2d-VEHvpa_2d,levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit3 = eBayes(fit2)
topTable(fit3,lfc=log2(1.5),p.value=0.05,adjust="BH",number=Inf)

