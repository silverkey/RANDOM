#!/usr/bin/R
library(RMySQL)
drv = dbDriver("MySQL")
db = 'mm9'
usr = 'mysql_dev'
pwd = 'dEvEl0pEr'
get.map = function(con) {
  query = paste(sep='',"SELECT DISTINCT name2,value ",
                       "FROM ensGene e, ensemblToGeneName g ",
                       "WHERE e.name=g.name")
  res = dbSendQuery(con,query)
  map = fetch(res,n=-1)
  colnames(map) = c('ens','sym')
  map
}
con = dbConnect(drv,db,usr,pwd)
map = get.map(con)
rm(drv,db,usr,pwd,get.map,con)
