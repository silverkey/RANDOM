#!/usr/local/bin/R

# Extract coord info, sample names and normalized signals from an OSC table
# The script is working on level2 OSC data table

# Naming conventions as used in the parvalbumin data
osc.name = 'level2.osc'
exp.field = '.H.'
norm.field = 'norm'
raw.field = 'raw'
del.1 = 'norm.H.'
del.2 = '_final'
del.3 = 'raw.H.'

# Read the OSC table
tab = read.delim(osc.name,comment.char="#")

# Column indexes to be used to split the table
exp.col = grep(exp.field,names(tab))
norm.col = grep(norm.field,names(tab)) 
raw.col = grep(raw.field,names(tab)) 

# Extract tables and give them proper rownames

# coord.table is the table with the coordinates info
coord.tab = tab[,-exp.col]

# norm.table is the table with the normalized expression values
norm.tab = tab[,norm.col]
colnames(norm.tab) = gsub(del.1,'',colnames(norm.tab))
colnames(norm.tab) = gsub(del.2,'',colnames(norm.tab))
norm.tab$id = coord.tab$id
order = c(ncol(norm.tab),seq(1,ncol(norm.tab)-1,1))
norm.tab = norm.tab[,order]

# raw.table is the table with the rawalized expression values
raw.tab = tab[,raw.col]
colnames(raw.tab) = gsub(del.3,'',colnames(raw.tab))
colnames(raw.tab) = gsub(del.2,'',colnames(raw.tab))
raw.tab$id = coord.tab$id
order = c(ncol(raw.tab),seq(1,ncol(raw.tab)-1,1))
raw.tab = raw.tab[,order]

# target is the table containing the info about the samples
target = data.frame(matrix(unlist(strsplit(colnames(norm.tab)[-1],'_')),ncol=3,byrow=T))
colnames(target) = c('chem','time','rep')
target$SampleName = colnames(norm.tab)[-1]
order = c(ncol(target),seq(1,ncol(target)-1,1))
target = target[,order]

# Write the 3 external tables to external txt files
write.table(coord.tab,file='coord.txt',sep="\t",quote=F,row.names=F)
write.table(norm.tab,file='norm.txt',sep="\t",quote=F,row.names=F)
write.table(raw.tab,file='raw.txt',sep="\t",quote=F,row.names=F)
write.table(target,file='target.txt',sep="\t",quote=F,row.names=F)
