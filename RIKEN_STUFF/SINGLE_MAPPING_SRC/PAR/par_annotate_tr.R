#!/usr/bin/R

library(RMySQL)

drv = dbDriver("MySQL")
db = 'mm9'
usr = 'mysql_dev'
pwd = 'dEvEl0pEr'

coord = read.delim(file='par_coord_filt.txt')
sel = read.delim(file='par_expdiff_results.xls')
out.file = 'tr_par_expdiff_results.xls'

get.overlap = function(chr,start,end) {
  query = paste(sep='',"SELECT DISTINCT name2,value FROM ensGene e, ensemblToGeneName g WHERE ",
                       "e.name=g.name AND ",
                       "chrom = '",chr,"' AND ",
                       "txStart <= ",end," AND ",
                       "txEnd >= ",start)
  res = dbSendQuery(con,query)
  pos = fetch(res,n=-1)
  pos
}

con = dbConnect(drv,db,usr,pwd)
common = unique(coord[unique(coord$id) %in% unique(sel$id),])
res = c()

for(i in 1:nrow(common)) {
  row = common[i,]
  id = as.character(row$id)
  chr = as.character(row$chr)
  start = as.character(row$start.0base)
  end = as.character(row$end)
  r = get.overlap(chr,start,end)
  if(nrow(r)>0) {
    for(ii in 1:nrow(r)) {
      a = c(id,as.character(r[ii,]))
      print(a)
    }
  }
  else {
    a = c(id,'NA','NA')
    print(a)
  }
  res = rbind(res,a)
}

rownames(res) = seq(1,nrow(res),1)
colnames(res) = c('id','ensgene','name')
res = as.data.frame(res)
annotated = merge(sel,res,all.x=T)
write.table(annotated,file=out.file,sep="\t",row.names=F,quote=F)
