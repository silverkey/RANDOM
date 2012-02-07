sel.dir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/send'
coord.file = '/Users/remo/Desktop/RIKEN/SINGLE_MAPPING/PYRAMIDAL/PYRAMIDAL/pyr_coord.txt'
id.col = 1
setwd(sel.dir)
tabs = dir()
data = tabs[grep('pr.+txt',tabs)]
names = gsub('.txt','.coord.txt',data)
coord = read.delim(file=coord.file)

for(i in 1:length(data)) {
  selected = read.delim(file=data[i])
  sel.coord = merge(selected,coord,by.x='ID',by.y='id',all.x='T')
  write.table(sel.coord,file=names[i],sep="\t",quote=F,row.names=F)
}

rm(list=ls())



sel.dir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/send'
coord.file = '/Users/remo/Desktop/RIKEN/SINGLE_MAPPING/PARVALBUMIN/PARVALBUMIN/par_coord.txt'
id.col = 1
setwd(sel.dir)
tabs = dir()
data = tabs[grep('pv.+txt',tabs)]
names = gsub('.txt','.coord.txt',data)
coord = read.delim(file=coord.file)

for(i in 1:length(data)) {
  selected = read.delim(file=data[i])
  sel.coord = merge(selected,coord,by.x='ID',by.y='id',all.x='T')
  write.table(sel.coord,file=names[i],sep="\t",quote=F,row.names=F)
}

rm(list=ls())
