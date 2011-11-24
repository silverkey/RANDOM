out.tab = 'repeat_enrichments_PYR.xls'

rep.name = as.character(rep$rep.family)
sig.rep.name = as.character(rep[rep$id %in% unique(TAB$id),]$rep.family)
unique.name = unique(as.character(rep[rep$id %in% unique(TAB$id),]$rep.family))
m1 = matrix(ncol=7,byrow=T)
colnames(m1) = c('rep','res.n','res.tot','pop.n','pop.tot','p.val','adj.p.val')
for(i in 1:length(unique.name)) {
  x = unique.name[i]
  a = sum(sig.rep.name==x)
  if(a >= 10) {
    b = sum(rep.name==x)
    r = prop.test(c(a,length(sig.rep.name)),c(b,length(rep.name)),alternative='g')
    mr = c(x,a,length(sig.rep.name),b,length(rep.name),r$p.value,NA)
    m1=rbind(m1,mr,deparse.level=0)
  }
}

rep.name = as.character(rep$rep.class)
sig.rep.name = as.character(rep[rep$id %in% unique(TAB$id),]$rep.class)
unique.name = unique(as.character(rep[rep$id %in% unique(TAB$id),]$rep.class))
m2 = matrix(ncol=7,byrow=T)
colnames(m2) = c('rep','res.n','res.tot','pop.n','pop.tot','p.val','adj.p.val')
for(i in 1:length(unique.name)) {
  x = unique.name[i]
  a = sum(sig.rep.name==x)
  if(a >= 10) {
    b = sum(rep.name==x)
    r = prop.test(c(a,length(sig.rep.name)),c(b,length(rep.name)),alternative='g')
    mr = c(x,a,length(sig.rep.name),b,length(rep.name),r$p.value,NA)
    m2=rbind(m2,mr,deparse.level=0)
  }
}

rep.name = as.character(rep$rep.name)
sig.rep.name = as.character(rep[rep$id %in% unique(TAB$id),]$rep.name)
unique.name = unique(as.character(rep[rep$id %in% unique(TAB$id),]$rep.name))
m3 = matrix(ncol=7,byrow=T)
colnames(m3) = c('rep','res.n','res.tot','pop.n','pop.tot','p.val','adj.p.val')
for(i in 1:length(unique.name)) {
  x = unique.name[i]
  a = sum(sig.rep.name==x)
  if(a >= 10) {
    b = sum(rep.name==x)
    r = prop.test(c(a,length(sig.rep.name)),c(b,length(rep.name)),alternative='g')
    mr = c(x,a,length(sig.rep.name),b,length(rep.name),r$p.value,NA)
    m3=rbind(m3,mr,deparse.level=0)
  }
}

m = rbind(m1[-1,],m2[-1,],m3[-1,])
m=as.data.frame(m)
m$adj.p.val = p.adjust(as.vector(m$p.val),method='BH')
sig.m = m[m$adj.p.val<=0.05,]
sig.m = sig.m[order(sig.m$adj.p.val),]

write.table(sig.m,file=out.tab,sep="\t",row.names=F,quote=F)

