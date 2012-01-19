#!/usr/bin/R

library(DESeq)
library(limma)
library(gplots)

out.file = 'pyr_expdiff_results.xls'
raw.file = 'pyr_raw_filt.txt'
tar.file = 'pyr_target.txt'
cl.file = 'pyr_DE_CLUSTERING.pdf'
cl.name = 'Pyramidal DE in'

min.reads = 25

raw = read.delim(file=raw.file)
raw = raw[rowSums(raw[,-1])>=min.reads,]
target = read.delim(file=tar.file)

# TSA_2d TSA_2h VEHtsa_2d VEHtsa_2h VEHvpa_2d VEHvpa_2h VPA_2d VPA_2h
TSA = c(1:12)
VPA = c(13:24)
ALL = c(TSA,VPA)

target$TSA_2h  = c(2,2,2,1,1,1,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0)
target$TSA_2d  = c(1,1,1,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0)
target$TSA_2hd = c(1,1,1,1,1,1,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0)

target$VPA_2h  = c(0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,1,1,1)
target$VPA_2d  = c(0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,1,1,1,2,2,2)
target$VPA_2hd = c(0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,1,1,1,1,1,1)

target$T.V_2h  = c(2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1)
target$T.V_2d  = c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2)
target$T.V_2hd = c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1)

target$chem_time = paste(target$chem,target$time,sep='_')

eset = raw[,-1]
rownames(eset) = raw$id

pdf(file=cl.file,paper='a4r',width=8.3,height=11.7,pointsize=8)

make.comparison = function(eset,group,name,chem) {
  
  eset = eset[,chem]
  group = group[group>0]

  cds = newCountDataSet(eset,group)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  vsd = getVarianceStabilizedData(cds)

  design = model.matrix(~0+factor(group))
  colnames(design) = gsub('factor(.+)1','A',colnames(design))
  colnames(design) = gsub('factor(.+)2','B',colnames(design))
  contrast = makeContrasts(A-B,levels=design)
  
  fit = lmFit(vsd,design=design)
  fit2 = contrasts.fit(fit,contrast)
  fit3 = eBayes(fit2)
  res = topTable(fit3,n='a',p.value=0.1,sort.by='logFC')  
  
  if(nrow(res) > 0) {
    if(nrow(res) > 1) {
      res.ordered = res[order(res$logFC,decreasing=T),]
      mat = as.matrix(vsd[res.ordered$ID,])
      heatmap.2(mat,col=redgreen(75),trace="none",dendrogram='none',
                Colv=F,Rowv=F,density.info='none',labRow=NA,margins=c(10,1),
                scale='row',colsep=seq(3,ncol(mat)-3,3),keysize=1,
                main=paste(cl.name,name,sep=' '))
    }
    res$comparison = name
    colnames(res) = gsub('ID','id',colnames(res))
    return(res[order(res$adj.P.Val),c('id','logFC','adj.P.Val','comparison')])
  }
  else return('NA')
}

name = 'TSA_2h'
group = target$TSA_2h
RES = make.comparison(eset,group,name,TSA)

name = 'TSA_2d'
group = target$TSA_2d
res = make.comparison(eset,group,name,TSA)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'TSA_2hd'
group = target$TSA_2hd
res = make.comparison(eset,group,name,TSA)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'VPA_2h'
group = target$VPA_2h
res = make.comparison(eset,group,name,VPA)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'VPA_2d'
group = target$VPA_2d
res = make.comparison(eset,group,name,VPA)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'VPA_2hd'
group = target$VPA_2hd
res = make.comparison(eset,group,name,VPA)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'T.V_2h'
group = target$T.V_2h
res = make.comparison(eset,group,name,ALL)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'T.V_2d'
group = target$T.V_2d
res = make.comparison(eset,group,name,ALL)
if(is.data.frame(res)) RES = rbind(RES,res)

name = 'T.V_2hd'
group = target$T.V_2hd
res = make.comparison(eset,group,name,ALL)
if(is.data.frame(res)) RES = rbind(RES,res)

dev.off()

write.table(RES,file=out.file,sep="\t",row.names=F,quote=F)
