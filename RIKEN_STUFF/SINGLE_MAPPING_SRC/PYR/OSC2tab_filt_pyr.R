#!/usr/local/bin/R

# Extract coord info, sample names and counts from an OSC table
# The script is working on level2 OSC data table and create also
# filtered tables

# Naming conventions as used in the pyramidal data
osc.name = 'pr_all_6_level2.osc'
exp.field = '_h_'
norm.field = 'norm'
raw.field = 'raw'
del.1 = 'norm.193_h_'
del.2 = '.notrim.+sorted'
del.3 = 'raw.193_h_'

# Cutoff on row sums
min.reads = 20

# Read the OSC table
tab = read.delim(osc.name,comment.char="#")

# Column indexes to be used to split the table
exp.col = grep(exp.field,names(tab),perl=T)
norm.col = grep(norm.field,names(tab),perl=T) 
raw.col = grep(raw.field,names(tab),perl=T) 

# Extract tables and give them proper rownames

# coord.table is the table with the coordinates info
coord.tab = tab[,-exp.col]

# norm.table is the table with the normalized expression values
norm.tab = tab[,norm.col]
colnames(norm.tab) = gsub(del.1,'',colnames(norm.tab),perl=T)
colnames(norm.tab) = gsub(del.2,'',colnames(norm.tab),perl=T)
norm.tab$id = coord.tab$id
order = c(ncol(norm.tab),seq(1,ncol(norm.tab)-1,1))
norm.tab = norm.tab[,order]

# raw.table is the table with the rawalized expression values
raw.tab = tab[,raw.col]
colnames(raw.tab) = gsub(del.3,'',colnames(raw.tab),perl=T)
colnames(raw.tab) = gsub(del.2,'',colnames(raw.tab),perl=T)
raw.tab$id = coord.tab$id
order = c(ncol(raw.tab),seq(1,ncol(raw.tab)-1,1))
raw.tab = raw.tab[,order]

# In the second chunck of OSC when we only use single mapped reads
# the colnames are different, here is a workaround to make them
# consistent with the past
colnames(norm.tab) = gsub('veh_','VEH',colnames(norm.tab),perl=T)
colnames(norm.tab) = gsub('^tsa','TSA',colnames(norm.tab),perl=T)
colnames(norm.tab) = gsub('^vpa','VPA',colnames(norm.tab),perl=T)

colnames(raw.tab) = gsub('veh_','VEH',colnames(raw.tab),perl=T)
colnames(raw.tab) = gsub('^tsa','TSA',colnames(raw.tab),perl=T)
colnames(raw.tab) = gsub('^vpa','VPA',colnames(raw.tab),perl=T)

# target is the table containing the info about the samples
target = data.frame(matrix(unlist(strsplit(colnames(norm.tab)[-1],'_')),ncol=3,byrow=T))
colnames(target) = c('chem','time','rep')
target$SampleName = colnames(norm.tab)[-1]
order = c(ncol(target),seq(1,ncol(target)-1,1))
target = target[,order]

# Write the 3 external tables to external txt files
#write.table(coord.tab,file='pyr_coord.txt',sep="\t",quote=F,row.names=F)
#write.table(norm.tab,file='pyr_norm.txt',sep="\t",quote=F,row.names=F)
#write.table(raw.tab,file='pyr_raw.txt',sep="\t",quote=F,row.names=F)
write.table(target,file='pyr_target.txt',sep="\t",quote=F,row.names=F)

# Filter based on total reads per row
rsum = rowSums(raw.tab[,-1])
write.table(coord.tab[rsum>=min.reads,],file='pyr_coord_filt.txt',sep="\t",quote=F,row.names=F)
#write.table(norm.tab[rsum>=min.reads,],file='pyr_norm_filt.txt',sep="\t",quote=F,row.names=F)
write.table(raw.tab[rsum>=min.reads,],file='pyr_raw_filt.txt',sep="\t",quote=F,row.names=F)

