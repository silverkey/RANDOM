# cond1 has to be the normal - the results will be refered to cond2 vs cond1
diff.exp = function(eset,target,target.comp.col,cond1,cond2) {

  snames = as.character(target[target.comp.col %in% c(cond1,cond2),'SampleName'])
  cnames = as.character(target.comp.col[target.comp.col %in% c(cond1,cond2)])  
  conds = factor(cnames)

  countsTable = eset[,colnames(eset) %in% snames]

  countsTable = countsTable[rowSums(countsTable)>100,]

  cds = newCountDataSet(countsTable,conds)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  res = nbinomTest(cds,cond1,cond2)
  nrow(res[res$padj<=0.1,])
	
  res[res$padj<=0.1,]
  countsTable[res[res$padj<=0.1,1],]

  selected = res[res$padj<=0.1,]
  selected = selected[order(selected$foldChange,decreasing=T),]
  selected$comparison = paste(cond2,'vs',cond1,sep='_')

  main = paste('DE in comparison',cond2,'vs',cond1)

  rank = c(10,11,12,4,5,6,7,8,9,1,2,3,16,17,18,22,23,24,13,14,15,19,20,21)

  heatmap.2(as.matrix(countsTable[selected$id,]),col=redgreen(75),trace="none",dendrogram='none',
            Colv=F,Rowv=F,density.info='none',labRow=NA,margins=c(10,1),
            scale='row',colsep=seq(3,ncol(raw[,-1])-3,3),keysize=1,main=main)

  heatmap.2(as.matrix(eset[selected$id,rank]),col=redgreen(75),trace="none",dendrogram='none',
            Colv=F,Rowv=F,density.info='none',labRow=NA,margins=c(10,1),
            scale='row',colsep=seq(3,ncol(raw[,-1])-3,3),keysize=1,main=main)

  selected
}

library(DESeq)
library(limma)
library(gplots)

rep.file = 'overlap_rna_rep_PAR.txt'
tr.file = 'simple_overlap_rna_transcripts_PAR.txt'
out.tab  = 'repeat_expdiff_results_PAR.xls'
expclfile = 'exp_clustering_PAR.pdf'
HCtitle = 'Hierarchical clustering dendrogram for raw counts in PARVALBUMIN'
cluster.file = 'CLUSTERING_PAR.pdf'
clname = 'PARVALBUMIN DE RE in'

# TSA_2d TSA_2h VEHtsa_2d VEHtsa_2h VEHvpa_2d VEHvpa_2h VPA_2d VPA_2h
TSA = c(4,5,6,10,11,12,1,2,3,7,8,9)
VPA = c(22,23,24,16,17,18,19,20,21,13,14,15)
ALL = c(TSA,VPA)

rep = read.delim(file=rep.file)
tr = read.delim(file=tr.file)
norm = read.delim(file='norm.txt')
raw = read.delim(file='raw.txt')
coord = read.delim(file='coord.txt')
eset = raw[,-1]
rownames(eset) = raw$id

target = read.delim(file='target.txt')
target$chem_time = paste(target$chem,target$time,sep='_')
target$summ_chem_time = c(rep('T_2d',3),rep('T_2h',3),rep('C_2d',3),rep('C_2h',3),rep('C_2d',3),rep('C_2h',3),rep('T_2d',3),rep('T_2h',3))

pdf(file=cluster.file,paper='a4r',width=8.3,height=11.7,pointsize=8)

target.comp.col = target$chem_time

cond1 = 'VEHtsa_2h'
cond2 = 'TSA_2h'
RES = diff.exp(eset,target,target.comp.col,cond1,cond2)

cond1 = 'VEHtsa_2d'
cond2 = 'TSA_2d'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

cond1 = 'VEHvpa_2h'
cond2 = 'VPA_2h'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

cond1 = 'VEHvpa_2d'
cond2 = 'VPA_2d'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

target.comp.col = target$summ_chem_time
cond1 = 'C_2h'
cond2 = 'T_2h'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

cond1 = 'C_2d'
cond2 = 'T_2d'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

target.comp.col = target$chem
cond1 = 'TSA'
cond2 = 'VPA'
res = diff.exp(eset,target,target.comp.col,cond1,cond2)
RES = rbind(RES,res)

dev.off()

RES.rep = merge(RES[,c(1,9,6,8)],rep[,c(1,5,14,15,16,17)],all.x=T)
RES.tr = merge(RES[,c(1,9,6,8)],tr[,c(1,5,8,7)],all.x=T)

write.table(RES,file=out.tab,sep="\t",row.names=F,quote=F)
write.table(RES.rep,file=paste('rep_',out.tab,sep=''),sep="\t",row.names=F,quote=F)
write.table(RES.tr,file=paste('tr_',out.tab,sep=''),sep="\t",row.names=F,quote=F)




raw[raw$id=='L2_level1_chr19_+_22011784',]
