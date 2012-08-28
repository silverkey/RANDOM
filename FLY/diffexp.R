x=read.table('Trinity_X.RSEM.fpkm',sep="\t",comment.char="#",colClasses="character",header=T,quote="")
y=read.table('Trinity_Y.RSEM.fpkm',sep="\t",comment.char="#",colClasses="character",header=T,quote="")

colnames(x)=c('tid','Xlength','Xeff_length','Xcount','Xfrac','Xfpkn','Xperc')
colnames(y)=c('tid','Ylength','Yeff_length','Ycount','Yfrac','Yfpkn','Yperc')

xy=merge(x,y,by.x='tid',by.y='tid')

write.table(xy,file='XY_RSEM_counts.tsv',sep="\t",row.names=F,quote=F)

interesting = na.omit(xy[as.numeric(xy$Ycount)>=100 & as.numeric(xy$Xcount)==0,c('tid','Xlength','Xeff_length','Xcount','Ycount','Xfpkn','Yfpkn')])

interesting = interesting[order(interesting$Ycount,decreasing=T),]

write.table(interesting,file='XY_RSEM_interesting.tsv',sep="\t",row.names=F,quote=F)
