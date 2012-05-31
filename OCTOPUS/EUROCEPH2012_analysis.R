/Users/remo/ANALYSIS/FIORITO/OCTOPUS/CAP3/FILTERED/IDMAP

c.i = read.delim('CAP2ID.txt')
m.c = read.delim('MIRA2CAP.txt')
m.r = read.delim('contigreadlist.parsed')

a1 = merge(m.r,m.c,by.x='id',by.y='MIRA',all.x=T)
a2 = merge(m.c,c.i,by.x='CAP3',by.y='CAP3',all.x=T)
a3 = merge(a1,a2,by.x='id',by.y='MIRA',all.x=T)
a4 = merge(a3,c.i,by.x='id',by.y='CAP3',all.x=T)

id = apply(a4[,7:8],1,function(x) sum(x,na.rm=T))
a4$final = id

cap = sapply(as.vector(a4$CAP3.x), function(x) ifelse(is.na(x),0,x),simplify = T)
a4$cap=as.character(cap)

res = a4[,c(1,10,9,2,3,4)]
colnames(res) = c('mira','cap','id','illumina','roche','sanger')

write.table(res,file='idmap.tab',sep="\t",row.names=F,quote=F)

#-----------
res = read.delim('idmap.tab')

anno = read.table('blast.anno.txt',sep="\t",comment.char="",colClasses="character",header=T,quote="")
anno$Uniref.Desc = gsub("\"","",anno$Uniref.Desc)
anno$CDD.Desc = gsub("\"","",anno$CDD.Desc)

tot = apply(res[,4:6],2,sum)

res$ill.prop = (res[,4]/as.numeric(tot['illumina']))*1000
res$roche.prop = (res[,5]/as.numeric(tot['roche']))*1000
res$sanger.prop = (res[,6]/as.numeric(tot['sanger']))*1000

#-----------
roche = res[res$roche>0,]

roche$roche.ill.p = apply(roche[,4:5],1,function(x) prop.test(c(x[1],x[2]),c(as.numeric(tot['illumina']),as.numeric(tot['roche'])))$p.value)
roche$roche.ill.pad = p.adjust(roche$roche.ill.p)
roche.ill.up = na.omit(roche[roche$roche.ill.pad<0.01 & roche$roche.prop/roche$ill.prop>=2 & roche$id>0,])
roche.ill.down = na.omit(roche[roche$roche.ill.pad<0.01 & roche$ill.prop/roche$roche.prop>=2 & roche$id>0,])

roche$roche.sanger.p = apply(roche[,5:6],1,function(x) prop.test(c(x[1],x[2]),c(as.numeric(tot['roche']),as.numeric(tot['sanger'])))$p.value)
roche$roche.sanger.pad = p.adjust(roche$roche.sanger.p)
roche.sanger.up = na.omit(roche[roche$roche.sanger.pad<0.01 & roche$roche.prop/roche$sanger.prop>=2 & roche$id>0,])
roche.sanger.down = na.omit(roche[roche$roche.sanger.pad<0.01 & roche$sanger.prop/roche$roche.prop>=2 & roche$id>0,])

write.table(roche.ill.up$id,file='roche.ill.up',sep="\t",row.names=F,quote=F)
write.table(unlist(strsplit(anno[anno$Identifier %in% roche.ill.up$id,'Uniref.Desc'],' n=.+$')),file='roche.ill.up.xls',sep="\t",row.names=F,quote=F)
write.table(roche.sanger.up$id,file='roche.sanger.up',sep="\t",row.names=F,quote=F)
write.table(unlist(strsplit(anno[anno$Identifier %in% roche.sanger.up$id,'Uniref.Desc'],' n=.+$')),file='roche.sanger.up.xls',sep="\t",row.names=F,quote=F)

#-----------
sanger = res[res$sanger>0,]

sanger$sanger.ill.p = apply(sanger[,4:6],1,function(x) prop.test(c(x[1],x[2]),c(as.numeric(tot['illumina']),as.numeric(tot['sanger'])))$p.value)
sanger$sanger.ill.pad = p.adjust(sanger$sanger.ill.p)
sanger.ill.up = na.omit(sanger[sanger$sanger.ill.pad<0.01 & sanger$sanger.prop/sanger$ill.prop>=2 & sanger$id>0,])
sanger.ill.down = na.omit(sanger[sanger$sanger.ill.pad<0.01 & sanger$ill.prop/sanger$sanger.prop>=2 & sanger$id>0,])

sanger$sanger.roche.p = apply(sanger[,5:6],1,function(x) prop.test(c(x[1],x[2]),c(as.numeric(tot['roche']),as.numeric(tot['sanger'])))$p.value)
sanger$sanger.roche.pad = p.adjust(sanger$sanger.roche.p)
sanger.roche.up = na.omit(sanger[sanger$sanger.roche.pad<0.01 & sanger$sanger.prop/sanger$roche.prop>=2 & sanger$id>0,])
sanger.roche.down = na.omit(sanger[sanger$sanger.roche.pad<0.01 & sanger$roche.prop/sanger$sanger.prop>=2 & sanger$id>0,])

write.table(sanger.ill.up$id,file='sanger.ill.up',sep="\t",row.names=F,quote=F)
write.table(unlist(strsplit(anno[anno$Identifier %in% sanger.ill.up$id,'Uniref.Desc'],' n=.+$')),file='sanger.ill.up.xls',sep="\t",row.names=F,quote=F)
write.table(sanger.roche.up$id,file='sanger.roche.up',sep="\t",row.names=F,quote=F)
write.table(unlist(strsplit(anno[anno$Identifier %in% sanger.roche.up$id,'Uniref.Desc'],' n=.+$')),file='sanger.roche.up.xls',sep="\t",row.names=F,quote=F)
#-----------
only.ill = res[res$ill.prop>=1 & res$roche.prop==0 & res$sanger==0 & res$id>0,]
only.roche = res[res$ill.prop==0 & res$roche.prop>=1 & res$sanger==0 & res$id>0,]
only.sanger = res[res$ill.prop==0 & res$roche.prop==0 & res$sanger>=1 & res$id>0,]

write.table(only.ill$id,file='only.ill',sep="\t",row.names=F,quote=F)
write.table(only.sanger$id,file='only.sanger',sep="\t",row.names=F,quote=F)
write.table(anno[anno$Identifier %in% only.ill$id,],file='only.ill.xls',sep="\t",row.names=F,quote=F)

#-----------
pdf(file='GO_composition.pdf',width=18,height=10)

go = read.table('blast2GO.ann')
names(go) = c('id','go')
godef = read.table('go_definition.txt',sep="\t",comment.char = "",colClasses = "character",quote="")
names(godef) = c('id','def','class')
goc = table(go$go)
top300 = names(sort(goc,decreasing=T)[1:300])

topF = godef[godef$id %in% top300 & godef$class=='F',]
tF = cbind(topF,goc[names(goc) %in% topF$id])
names(tF) = c('id','def','class','n')
tF = na.omit(tF[order(tF$n,decreasing=T)[1:15],])
pie(tF$n,labels=tF$def,cex=1.5,cex.main=2,main='Top 15 GO Molecular Function')
tF = tF[order(tF$n),]
dotchart(tF$n,labels=tF$def,cex=1.8,main='Top 15 GO Molecular Function')

topC = godef[godef$id %in% top300 & godef$class=='C',]
tC = cbind(topC,goc[names(goc) %in% topC$id])
names(tC) = c('id','def','class','n')
tC = na.omit(tC[order(tC$n,decreasing=T)[1:15],])
pie(tC$n,labels=tC$def,cex=1.5,cex.main=2,main='Top 15 GO Cellular Component')
tC = tC[order(tC$n),]
dotchart(tC$n,labels=tC$def,cex=1.8,main='Top 15 GO Cellular Component')

topP = godef[godef$id %in% top300 & godef$class=='P',]
tP = cbind(topP,goc[names(goc) %in% topP$id])
names(tP) = c('id','def','class','n')
tP = na.omit(tP[order(tP$n,decreasing=T)[1:15],])
pie(tP$n,labels=tP$def,cex=1.5,cex.main=2,main='Top 15 GO Biological Process')
tP = tP[order(tP$n),]
dotchart(tP$n,labels=tP$def,cex=1.8,main='Top 15 GO Biological Process')

dev.off()

#anno[anno$Identifier %in% unique(go[go$go=='GO:0002119',1]),]$Uniref.Desc

#-----------
pdf(file='domain.pdf',width=18,height=10)

cdd = table(anno$CDD.Desc)
top100 = names(sort(cdd,decreasing=T)[1:100])
topCDD = cdd[names(cdd) %in% top100[1:10]]

par(las=2) # make label text perpendicular to axis
par(mar=c(5,22,4,2)) # increase y-axis margin.
dn = c('immunoglobulin domains','serpentine chemoreceptor','papillomavirus E5','collagen','transposase','reverse transcriptase','ankyrin repeats','WD40 domain','Zn-finger')
#barplot(sort(topCDD)[1:9],horiz=T,names=dn,cex.names=1.5,cex.axis=1.5)
barplot(sort(topCDD)[5:9],horiz=T,names=dn[5:9],cex.names=2.5,cex.axis=2,main='Top 5 Protein Domains',cex.main=2)

dev.off()


#-----------
pdf(file='noncoding.pdf',width=18,height=10)

res = res[res$id>0,]

tot = apply(res[,4:6],2,sum)

res$ill.prop = (res[,4]/as.numeric(tot['illumina']))*1000
res$roche.prop = (res[,5]/as.numeric(tot['roche']))*1000
res$sanger.prop = (res[,6]/as.numeric(tot['sanger']))*1000

noncoding = res[res$id %in% anno[anno$Coding==0,]$Identifier,]
noncoding = res[res$id %in% anno[anno$Coding==0 & anno$Portrait.Score>0.5,]$Identifier,]
nrow(anno[anno$Coding==0 & anno$Portrait.Score>0.5,])
nc.illumina = nrow(noncoding[noncoding$illumina>=1,])/tot['illumina']*100
nc.roche = nrow(noncoding[noncoding$roche>=1,])/tot['roche']*100
nc.sanger = nrow(noncoding[noncoding$sanger>=1,])/tot['sanger']*100

par(las=2) # make label text perpendicular to axis
par(mar=c(5,12,4,2)) # increase y-axis margin.
barplot(c(nc.illumina,nc.roche,nc.sanger),names=c('organism','neural','brain'),horiz=T,cex.names=3,cex.axis=2,main='Percentage Noncoding',cex.main=2)
dev.off()



#-----------
pdf(file='repeats.pdf',width=18,height=10)

rep = read.table('REP.tab',sep="\t")
names(rep) = c('id','rep','desc','start','end')
id.rep = table(rep[,1:2])
rname = colnames(id.rep)
repres = c()
for(i in 1:length(rname)) {
  n = rname[i]
  r = id.rep[,i]
  pos = length(r[r>0])
  repres[i] = pos
}
tr = c() #data.frame()
tr$rep = rname
tr$n = repres
tr = as.data.frame(tr)

transposon =id.rep[,c(8:12,16)]
transposon.pos = nrow(transposon[rowSums(transposon)>0,])

sine = id.rep[,13:15]
sine.pos = nrow(sine[rowSums(sine)>0,])

lsu = id.rep[,6:7]
lsu.pos = nrow(lsu[rowSums(lsu)>0,])

other = id.rep[,c(1:5)]
other.pos = nrow(other[rowSums(other)>0,])

par(las=2) # make label text perpendicular to axis
par(mar=c(5,15,4,2)) # increase y-axis margin.
barplot(c(other.pos,lsu.pos,transposon.pos,sine.pos),
        horiz=T,names=c('other','lsu','transposon','SINE'),main='Number of Transcripts Matching Repeats',
        cex.names=3,cex.axis=2,cex.main=2)

dev.off()














