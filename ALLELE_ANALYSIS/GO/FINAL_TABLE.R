fbp = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/BP/alleles_GOBP_significant.xls'
fmf = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/MF/alleles_GOMF_significant.xls'
fcc = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/CC/alleles_GOCC_significant.xls'

get.tab = function(fname) {
  tab = read.table(file=fname,sep='\t',header=T)
  tab
}

bp = get.tab(fbp)
bp$division = 'BP'

mf = get.tab(fmf)
mf$division = 'MF'

#cc = get.tab(cc)
#cc$division = 'cc'

sel = rbind(bp[bp$adj_p<=0.1,],mf[mf$adj_p<=0.1,])

ratio = ((sel$All_in_class/sel$All_in_GO)*100)/((sel$No_all_in_class/sel$No_all_in_GO)*100)

for(i in 1:length(ratio)) {
  ratio[i] = ifelse(ratio[i] > 1, ratio[i], -1/ratio[i])
}

sel$ratio = ratio

write.table(sel,file='FC_alleles_GO_significant.xls',sep='\t',row.names=F,quote=F)


fbp = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/BP/alleles_GOBP_significant.xls'
fmf = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/MF/alleles_GOMF_significant.xls'
fcc = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/GO_ALLELES/NEW/CC/alleles_GOCC_significant.xls'

get.tab = function(fname) {
  tab = read.table(file=fname,sep='\t',header=T)
  tab
}

bp = get.tab(fbp)
bp$division = 'BP'

mf = get.tab(fmf)
mf$division = 'MF'

#cc = get.tab(cc)
#cc$division = 'cc'

sel = rbind(bp,mf)

ratio = ((sel$All_in_class/sel$All_in_GO)*100)/((sel$No_all_in_class/sel$No_all_in_GO)*100)

for(i in 1:length(ratio)) {
  ratio[i] = ifelse(ratio[i] > 1, ratio[i], -1/ratio[i])
}

sel$ratio = ratio

write.table(sel,file='FC_alleles_GO_all.xls',sep='\t',row.names=F,quote=F)

