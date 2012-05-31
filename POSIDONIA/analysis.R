#!/usr/bin/R

# Script by Remo Sanges 21/05/2012

# reading the files
a = read.table('pooc_a_go.txt',sep='\t',comment.char="",colClasses="character",header=T,quote="")
b = read.table('pooc_b_go.txt',sep='\t',comment.char="",colClasses="character",header=T,quote="")
cloneid = read.table('array_id.txt',sep='\t',comment.char="",colClasses="character",header=T,quote="")
selected = read.table('array_selected.txt',sep='\t',comment.char="",colClasses="character",header=T,quote="")

# adjusting the ids
b$cloneid = gsub('Contig','Pooc_B_c',b$cloneid)

# concatenate A + B
c = rbind(a[,c(1,3)],b)

# solving repetition in one file
cloneid = unique(cloneid)
cloneid = as.data.frame(cloneid[1:977,])
colnames(cloneid) = 'cloneid'

# make the universe
clonego = merge(cloneid,c,by.x='cloneid',by.y='cloneid',all.x=T,all.y=T)
clonego = na.omit(clonego)

# selected ids
all = as.data.frame(selected$Name)
up = as.data.frame(selected[selected$regulation == 'up','Name'])
down = as.data.frame(selected[selected$regulation == 'down','Name'])

# write tables
write.table(all,file='all.txt',sep="\t",row.names=F,quote=F)
write.table(up,file='up.txt',sep="\t",row.names=F,quote=F)
write.table(down,file='down.txt',sep="\t",row.names=F,quote=F)
write.table(clonego,file='GO_array_annotation.txt',sep="\t",row.names=F,quote=F)
