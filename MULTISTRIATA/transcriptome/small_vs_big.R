library(limma)
library(genefilter)
target = readTargets('CLASSES.xls')
target$mat.size = paste(target$Mating,target$Size,sep='.')

brh = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/BRH_30.txt')
#brh = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/BRH/PM_SR_tblastx/blast1.table')
#brh = brh[brh$coverage>=50 & brh$evalue<=0.0000000001,1:2]
#colnames(brh) = c('PM','SR')

#pm = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/multistriata/RNA/Pseudo-nitzschia_multistriata_counting/counts.txt')
#sr = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/robusta/RNA/counting/counts.txt')
pm = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/multistriata/RNA/Pseudo-nitzschia_multistriata_counting/normal_counts.txt')
sr = read.delim('/Users/remo/ANALYSIS/multistriata/JGI/robusta/RNA/counting/normal_counts.txt')

anno = read.table('/Users/remo/ANALYSIS/multistriata/JGI/PM_annotations.xls',sep="\t",comment.char="",colClasses="character",header=T,quote="")
pm.anno = merge(pm,anno,by.x='Gene',by.y='Identifier',all.x=T,sort=F)

pm.anno.filt = pm.anno[pm.anno$Rfam.Desc=='NULL',]
pm.anno.filt = pm.anno.filt[grep('ribosomal protein',pm.anno.filt$Uniref.Desc,invert=T,fixed=T),]
pm.anno.filt = pm.anno.filt[grep('Ribosomal protein',pm.anno.filt$Uniref.Desc,invert=T,fixed=T),]
pm.anno.filt = pm.anno.filt[grep('ribosomal protein',pm.anno.filt$CDD.Desc,invert=T,fixed=T),]
pm.anno.filt = pm.anno.filt[grep('Ribosomal protein',pm.anno.filt$CDD.Desc,invert=T,fixed=T),]

pm.filt = pm[pm$Gene %in% pm.anno.filt$Gene,]
brh.filt = brh[brh$PM %in% pm.anno.filt$Gene,]

a = merge(brh.filt,pm.filt,by.x='PM',by.y='Gene',sort=F)[,c(1,2,3,5,6)]
b = merge(a,sr,by.x='SR',by.y='Gene',sort=F)[,c(1,2,3,4,5,6,8,9,10,11)]
c = b
names(c) = c('SR','PM','pm.length','pm.CIIO','pm.CIIP','sr.length','sr.CIIG','sr.CIIH','sr.CIII','sr.CIIN')
t = c[,c(1,2,3,6,4,5,7,8,9,10)]
cor(t[,5:10])

exprs = t[,c(2,5:10)]
names(exprs) = c('PM','CIIO','CIIP','CIIG','CIIH','CIII','CIIN')

#-----------------------
# BIG - SMALL ANALYSIS
#-----------------------
design = model.matrix(~0+target$Size)
colnames(design) = gsub('target$Size','',colnames(design),fixed=T)
contrast.matrix = makeContrasts(small-big,levels=design)

eset = exprs[,2:7]
rownames(eset) = exprs[,1]

f1 = pOverA(0.5,1)
f2 = function(x) (IQR(x) > 0.5)
ff = filterfun(f1, f2)
selected = genefilter(eset, ff)
eset = eset[selected,]

fit = lmFit(eset,design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit3 = eBayes(fit2)
res = topTable(fit3,n=15)[,c(1,2,5,6)]
res = merge(res,anno[,c(2,5,9,13,17,18,20,22,24)],by.x='ID',by.y='Identifier',all.x=T,sort=F)
res = merge(res,brh,by.x='ID',by.y='PM',all.x=T,sort=F)
write.table(res,file='small_vs_big_top15.xls',sep="\t",row.names=F,quote=F)

#pdf(width=8,height=10,pointsize=1,bg="transparent",file='plus_vs_minus_top15.pdf')
pdf(bg="transparent",width=8,height=10,pointsize=1,paper='a4',file='small_vs_big_top15.pdf')
par(mfrow=c(4,2))
for(i in 1:nrow(res)) {
  ID=res[i,1]
  ann = anno[anno$Identifier==ID,]
  main = c(substr(ann[,'Uniref.Desc'],1,80),substr(ann[,'CDD.Desc'],1,80))
  dotchart(as.numeric(eset[ID,]),labels=colnames(eset),color=c(2,2,1,2,1,2),main=main,cex.main=0.9,cex=1.05)
}
dev.off()
