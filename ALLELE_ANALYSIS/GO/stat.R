pcut = 0.1
ncut = 100

o = read.table(file='alleles_GOBP.txt',sep='\t')
t = read.table(file='alleles_GOBP.txt',sep='\t')
colnames(t) = c('Class','All_in_GO','No_all_in_GO','All_in_class','No_all_in_class')
t = t[(t[,4]+t[,5])>=ncut,]

a = c()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[4],r[5])),as.numeric(c(r[2],r[3])))$p.value
  a = c(a,pval)
}

t$p_val = a
t$adj_p = p.adjust(a,'fdr')
t[t$adj_p<=pcut,]

write.table(t,file='alleles_GOBP_significant.xls',sep='\t',row.names=F,quote=F)

