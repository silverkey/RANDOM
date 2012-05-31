#!/usr/bin/R
library(RMySQL)

drv = dbDriver("MySQL")
db = 'mm9'
usr = 'genome'
pwd = ''
host = 'genome-mysql.cse.ucsc.edu'

get.sym.map = function(con) {
  query = paste(sep='',"SELECT DISTINCT name2,value ",
                "FROM ensGene e, ensemblToGeneName g ",
                "WHERE e.name=g.name")
  res = dbSendQuery(con,query)
  map = fetch(res,n=-1)
  colnames(map) = c('ens','sym')
  map
}

get.tra.map = function(con) {
  query = paste(sep='',"SELECT DISTINCT gene,transcript FROM ensGtp")
  res = dbSendQuery(con,query)
  map = fetch(res,n=-1)
  colnames(map) = c('ens','tra')
  map
}

get.bio.map = function(con) {
  query = paste(sep='',"SELECT DISTINCT name,source FROM ensemblSource")
  res = dbSendQuery(con,query)
  map = fetch(res,n=-1)
  colnames(map) = c('tra','bio')
  map
}

con = dbConnect(drv,db,usr,pwd,host)
sym.map = get.sym.map(con)
tra.map = get.tra.map(con)
bio.map = get.bio.map(con)

rm(drv,db,usr,pwd,host,get.sym.map,get.tra.map,get.bio.map,con)
