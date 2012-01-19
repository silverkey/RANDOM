#!/usr/bin/R

library(RMySQL)
library(hash)

drv = dbDriver("MySQL")
db = 'mm9'
usr = 'mysql_dev'
pwd = 'dEvEl0pEr'
in.file = 'tr_rep_par_expdiff_results.xls'
out.file = 'repeat_enrichments_PAR.xls'
min.n = 10

# Collect all the names of the tables containing the mapping of repeats
get.rep.tabs = function(con) {
  query = 'SHOW TABLES'
  res = dbSendQuery(con,query)
  res = fetch(res,n=-1)
  tabs = res[grep('rmsk',res[,1]),1]
  tabs
}

# Prepare the data with the total number of repeats in the genome
count.tot.rep = function(con,tabs,field){
  tot = hash()
  for(i in 1:length(tabs)) {
    query = paste(sep='',"SELECT DISTINCT ",field,", COUNT(*) ",
                  "FROM ",tabs[i], " GROUP BY ",field)
    res = dbSendQuery(con,query)
    res = fetch(res,n=-1)
    for(ii in 1:nrow(res)) {
      if(has.key(res[ii,1],tot)) {
        tot[[res[ii,1]]] = tot[[res[ii,1]]] + res[ii,2]
      }
      else {
        tot[[res[ii,1]]] = res[ii,2]
      }
    }
  }
  tot
}

# Test the enrichment of each specific kind of repeat
test.enrichment = function(res.rep,gen.rep,field) {
  m = matrix(ncol=8,byrow=T)
  colnames(m) = c('rep','res.n','res.tot','pop.n','pop.tot','p.val','adj.p.val','field')
  names = unique(res.rep)
  for(i in 1:length(names)) {
    res.n = sum(res.rep == names[i])
    res.tot = length(res.rep)
    pop.n = gen.rep[[names[i]]]
    pop.tot = sum(values(gen.rep))
    p.val = prop.test(c(res.n,pop.n),c(res.tot,pop.tot),alternative='g')$p.value
    row = c(names[i],res.n,res.tot,pop.n,pop.tot,p.val,NA,field)
    m = rbind(m,row,deparse.level=0)
  }
  m
}

con = dbConnect(drv,db,usr,pwd)
rep.tabs = get.rep.tabs(con)

tot.name = count.tot.rep(con,rep.tabs,'repName')
tot.class = count.tot.rep(con,rep.tabs,'repClass')
tot.family = count.tot.rep(con,rep.tabs,'repFamily')

res.tab = read.delim(file=in.file)
rep.name = as.vector(na.omit(res.tab$repName))
rep.class = as.vector(na.omit(res.tab$repClass))
rep.family = as.vector(na.omit(res.tab$repFamily))

m.name = test.enrichment(rep.name,tot.name,'repName')
m.class = test.enrichment(rep.class,tot.class,'repClass')
m.family = test.enrichment(rep.family,tot.family,'repFamily')

m = rbind(m.name[-1,],m.class[-1,],m.family[-1,])
m = as.data.frame(m)
m$adj.p.val = p.adjust(as.vector(m$p.val),method='BH')
sig.m = m[m$adj.p.val<=0.05,]
sig.m = sig.m[order(sig.m$adj.p.val),]

filtered = sig.m[as.numeric(as.character(sig.m$res.n))>=min.n,]
write.table(filtered,file=out.file,sep="\t",row.names=F,quote=F)
