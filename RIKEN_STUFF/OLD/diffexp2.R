#!/usr/bin/R

rep.file = 'overlap_rna_rep_PYR.txt'
out.tab  = 'repeat_expdiff_results_PYR.txt'

library(limma)

make.comparison =function(contrast,fit) {
  contrast.matrix = makeContrasts(contrast,levels=design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit3 = eBayes(fit2)
  res = topTable(fit3,adjust="BH",number=Inf)
  res$comparison = contrast
  res
}

plot.region = function(id,colnames) {
  data = matrix(as.numeric(exp.rep[exp.rep$id==id,][-1]),ncol=length(colnames))
  colnames(data) = colnames
  boxplot(data,col=8,cex.axis=0.6,main=paste(id,'normalized counts'))
}

rep = read.delim(file=rep.file)
norm = read.delim(file='norm.txt')
raw = read.delim(file='raw.txt')
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

# CHANGE AGAINST CONTROL
contrast = 'TSA_2h-VEHtsa_2h'
RES = make.comparison(contrast,fit)

contrast = 'TSA_2d-VEHtsa_2d'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

contrast = 'VPA_2h-VEHvpa_2h'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

contrast = 'VPA_2d-VEHvpa_2d'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

# CHANGE BETWEEN TREATMENT
contrast = '(TSA_2h-VEHtsa_2h)-(VPA_2h-VEHvpa_2h)'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

contrast = '(TSA_2d-VEHtsa_2d)-(VPA_2d-VEHvpa_2d)'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

# CHANGE BETWEEN TIME
contrast = '(TSA_2d-VEHtsa_2d)-(TSA_2h-VEHtsa_2h)'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

contrast = '(VPA_2d-VEHvpa_2d)-(VPA_2h-VEHvpa_2h)'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

contrast = '((TSA_2d-VEHtsa_2d)-(TSA_2h-VEHtsa_2h))-((VPA_2d-VEHvpa_2d)-(VPA_2h-VEHvpa_2h))'
res = make.comparison(contrast,fit)
RES = rbind(RES,res)

# CHANGE BETWEEN VEHICLE
#contrast = 'VEHtsa_2h-VEHvpa_2h'
#res = make.comparison(contrast,fit)
#RES = rbind(RES,res)
#
#contrast = 'VEHtsa_2d-VEHvpa_2d'
#res = make.comparison(contrast,fit)
#RES = rbind(RES,res)

# SIGNIFICANT
SIG = RES[RES$adj.P.Val<=0.05,]
TAB = merge(SIG[,c(1,2,6,8)],rep[,c(1,2,3,4,5,8,14,15,16,17)],by.x='ID', by.y='id', all.x=T)
TAB = merge(TAB,norm,by.x='ID',by.y='id',all.x=T)
write.table(TAB,file=out.tab,sep="\t",row.names=F,quote=F)
