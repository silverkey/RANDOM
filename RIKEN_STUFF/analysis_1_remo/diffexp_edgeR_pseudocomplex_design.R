#!/usr/bin/R

library(edgeR)
library(limma)
library(gplots)

rep.file = 'overlap_rna_rep_PYR.txt'
out.tab  = 'repeat_expdiff_results_PYR.xls'
expclfile = 'exp_clustering_PYR.pdf'
HCtitle = 'Hierarchical clustering dendrogram for raw counts in PYRAMIDAL'
MDStitle = 'MDS clustering of PYRAMIDAL samples'
cluster.file = 'CLUSTERING_PYR.pdf'
clname = 'PYRAMIDAL DE RE in'

TSA = c(4,5,6,10,11,12,1,2,3,7,8,9)
VPA = c(22,23,24,16,17,18,19,20,21,13,14,15)
ALL = c(TSA,VPA)

rep = read.delim(file=rep.file)
norm = read.delim(file='norm.txt')
raw = read.delim(file='raw.txt')
coord = read.delim(file='coord.txt')
exp.rep = raw
eset = exp.rep[,-1]
rownames(eset) = exp.rep$id

# TSA_2d TSA_2h VEHtsa_2d VEHtsa_2h VEHvpa_2d VEHvpa_2h VPA_2d VPA_2h

target = read.delim(file='target.txt')

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

# Exp clustering
pdf(file=expclfile,paper='a4r',width=8.3,height=11.7,pointsize=8,title=MDStitle)
ecTr = dist(t(eset), method = "euclidean")
hecTr = hclust(ecTr, method = "average")
plot(hecTr, main = HCtitle, xlab = "")
# MDS generation
dd = DGEList(eset,group=target$chem_time)
col = c(2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9)
plotMDS.dge(dd,col=col,main=MDStitle)
dev.off()
  
pdf(file=cluster.file,paper='a4r',width=8.3,height=11.7,pointsize=8)

make.comparison = function(eset,group,rep,name,chem) {

  d = DGEList(eset,group=as.character(group))
  d = d[rowSums(1e+06 * d$counts/expandAsMatrix(d$samples$lib.size, dim(d)) > 10) >= 3, ]
  d = calcNormFactors(d)
  d = estimateCommonDisp(d)

  # It seems the calculation is done so the FC = 1-2
  # Code checked it does the second minus the first
  et = exactTest(d,pair=as.character(c('2','1')))

  stat = topTags(et,n=Inf)
  tab = stat$table
  sig.tab = tab[tab$adj.P.Val<=0.05,]
  res = sig.tab[rownames(sig.tab) %in% rep$id,]
  res$comparison = name
  
  res.ordered = res[order(res$logFC,decreasing=T),]
  fitted = d$pseudo.alt
  mat = as.matrix(fitted[rownames(res.ordered),])
  heatmap.2(mat[,chem],col=redgreen(75),trace="none",dendrogram='none',
            Colv=F,Rowv=F,density.info='none',labRow=NA,margins=c(10,1),
            scale='row',colsep=seq(3,ncol(mat)-3,3),keysize=1,
            main=paste(clname,name,sep=' '))

  res$id = rownames(res)
  res[order(res$adj.P.Val),]
}

name = 'TSA_2h'
group = target$TSA_2h
RES = make.comparison(eset,group,rep,name,TSA)

name = 'TSA_2d'
group = target$TSA_2d
res = make.comparison(eset,group,rep,name,TSA)
RES = rbind(RES,res)

name = 'TSA_2hd'
group = target$TSA_2hd
res = make.comparison(eset,group,rep,name,TSA)
RES = rbind(RES,res)

name = 'VPA_2h'
group = target$VPA_2h
res = make.comparison(eset,group,rep,name,VPA)
RES = rbind(RES,res)

name = 'VPA_2d'
group = target$VPA_2d
res = make.comparison(eset,group,rep,name,VPA)
RES = rbind(RES,res)

name = 'VPA_2hd'
group = target$VPA_2hd
res = make.comparison(eset,group,rep,name,VPA)
RES = rbind(RES,res)

name = 'T.V_2h'
group = target$T.V_2h
res = make.comparison(eset,group,rep,name,ALL)
RES = rbind(RES,res)

name = 'T.V_2d'
group = target$T.V_2d
res = make.comparison(eset,group,rep,name,ALL)
RES = rbind(RES,res)

name = 'T.V_2hd'
group = target$T.V_2hd
res = make.comparison(eset,group,rep,name,ALL)
RES = rbind(RES,res)

dev.off()

# SIGNIFICANT
SIG = RES[RES$adj.P.Val<=0.05,]
TAB = merge(SIG[,c(6,2,4,5)],rep[,c(1,2,3,4,5,8,14,15,16,17)],by.x='id', by.y='id', all.x=T)
write.table(TAB,file=out.tab,sep="\t",row.names=F,quote=F)
