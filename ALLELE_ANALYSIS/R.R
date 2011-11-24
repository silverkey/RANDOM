pdf(file='RES.pdf',paper='a4r',width=8.3,height=11.7,pointsize=5)
par(mfrow=c(2,1))

# -----------
# BOXPLOT AND XYPLOT IN EACH POSITION
res = list()
tab = read.delim(file='alleles_scan_for_R.tab',header=F)
for(i in unique(tab[,3])) {
  cat(i)
  c = tab[tab[,3]==i,4]
  res[[as.character(i)]] = c
}

boxplot(res,outline=F,col=8,
        main='Percentage identity in specific position intervals',
        ylab='percentage identity',xlab='position')

plot(seq(-490,500,10),lapply(res,mean),type='b',ylim=c(55,86),
     main='Average identity between alleles genomic regions',
     xlab='position',ylab='average percentage identity')

# BOXPLOT AND T.TEST IN 2 POSITIONS (PRE, POST TSS) -100/100
# -----------
tab = read.table(file='FC_alleles_analysis_CGI2/tss100.tab',sep="\t")
boxplot(tab[,2],tab[,3],ylim=c(0,100),col=8,horizontal=T,outline=F,
        names=c('before TSS','after TSS'),
        main='Percentage identity between alleles genomic regions -100/100',
        xlab='percentage identity')

t.test(tab[,2],tab[,3])
t.test(tab[,2],tab[,3])$p.value

# BOXPLOT AND T.TEST IN 2 POSITIONS (PRE, POST TSS) -500/500
# -----------
tab = read.table(file='FC_alleles_analysis_CGI2/tss500.tab',sep="\t")

boxplot(tab[,2],tab[,3],ylim=c(0,100),col=8,horizontal=T,outline=F,
        names=c('before TSS','after TSS'),
        main='Percentage identity between alleles genomic regions -500/500',
        xlab='percentage identity')

t.test(tab[,2],tab[,3])
t.test(tab[,2],tab[,3])$p.value

# TRIPLETS CHANGES
# -----------
w = na.omit(read.delim(file='FC_alleles_analysis_CGI3/word100.tab',header=F))
w = w[w[, 2] != "NNN", ]
names(w) = c('zone','word','conserved','notconserved')
w$ratio = w$notconserved/w$conserved

pre.w = w[w$zone=='pre',]
pre.w$m = mean(pre.w$ratio)
pre.w$s = sd(pre.w$ratio)
pre.w$Z = (pre.w$ratio - pre.w$m) / pre.w$s
pre.w$pval = 1-pnorm(abs(pre.w$Z))
pre.w$adjp = p.adjust(pre.w$pval,method='BH')
pre.w[abs(pre.w$Z)>=2,]

post.w = w[w$zone=='post',]
post.w$m = mean(post.w$ratio)
post.w$s = sd(post.w$ratio)
post.w$Z = (post.w$ratio - post.w$m) / post.w$s
post.w$pval = 1-pnorm(abs(post.w$Z))
post.w$adjp = p.adjust(post.w$pval,method='BH')
post.w[abs(post.w$Z)>=2,]

par(cex=0.8,lwd=0.1)
plot(seq(1:nrow(pre.w)),pre.w$ratio,type='n',axes=F,ylim=c(0.55,0.91),
     xlab='',ylab='not-conserved / conserved ratio',main='Trinucleotide changes -100 / 0 bp before TSS')
rect(0.3,pre.w$m-(2*pre.w$s),nrow(pre.w)+2,pre.w$m+(2*pre.w$s),col='lightgray',border=NA)
rect(0.3,0.55,nrow(pre.w)+2,pre.w$m-(2*pre.w$s),col='lightpink',border=NA)
rect(0.3,pre.w$m+(2*pre.w$s),nrow(pre.w)+2,0.91,col='lightgreen',border=NA)
axis(2)
text(seq(1:nrow(pre.w)),pre.w$ratio,pre.w$word)
abline(h=pre.w$m-(2*pre.w$s),col='blue')
abline(h=pre.w$m+(2*pre.w$s),col='blue')
abline(h=pre.w$m,col='red')
#abline(h=0.55)
#abline(h=0.91)

par(cex=0.8,lwd=0.1)
plot(seq(1:nrow(post.w)),post.w$ratio,type='n',axes=F,ylim=c(0.45,0.67),
     xlab='',ylab='not-conserved / conserved ratio',main='Trinucleotide changes 0 / +100 bp after TSS')
rect(0.3,post.w$m-(2*post.w$s),nrow(post.w)+2,post.w$m+(2*post.w$s),col='lightgray',border=NA)
rect(0.3,0.45,nrow(post.w)+2,post.w$m-(2*post.w$s),col='lightpink',border=NA)
rect(0.3,post.w$m+(2*post.w$s),nrow(post.w)+2,0.67,col='lightgreen',border=NA)
axis(2)
text(seq(1:nrow(post.w)),post.w$ratio,post.w$word)
abline(h=post.w$m-(2*post.w$s),col='blue')
abline(h=post.w$m+(2*post.w$s),col='blue')
abline(h=post.w$m,col='red')
#abline(h=0.45)
#abline(h=0.67)

# -----------
dev.off()
