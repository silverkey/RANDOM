ffam = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/PFAM/FAMILY/alleles_PFAM_FAMILY_significant.xls'

fdom = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/PFAM/DOMAIN/alleles_PFAM_DOMAIN_significant.xls'

get.tab = function(fname) {
  tab = read.table(file=fname,sep='\t',header=T)
  tab
}

fam = get.tab(ffam)
fam$division = 'family'

dom = get.tab(fdom)
dom$division = 'domain'

sel = rbind(fam[fam$adj_p<=0.1,],dom[dom$adj_p<=0.1,])

ratio = ((sel$All_in_class/sel$All_in_PFAM)*100)/((sel$No_all_in_class/sel$No_all_in_PFAM)*100)

for(i in 1:length(ratio)) {
  ratio[i] = ifelse(ratio[i] > 1, ratio[i], -1/ratio[i])
}

sel$ratio = ratio

write.table(sel,file='FC_alleles_PFAM_significant.xls',sep='\t',row.names=F,quote=F)




ffam = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/PFAM/FAMILY/alleles_PFAM_FAMILY_significant.xls'

fdom = '/Users/remo/ANALYSIS/MOCK/FC_ALLELES/ANALYSIS_2/PFAM/DOMAIN/alleles_PFAM_DOMAIN_significant.xls'

get.tab = function(fname) {
  tab = read.table(file=fname,sep='\t',header=T)
  tab
}

fam = get.tab(ffam)
fam$division = 'family'

dom = get.tab(fdom)
dom$division = 'domain'

sel = rbind(fam,dom)

ratio = ((sel$All_in_class/sel$All_in_PFAM)*100)/((sel$No_all_in_class/sel$No_all_in_PFAM)*100)

for(i in 1:length(ratio)) {
  ratio[i] = ifelse(ratio[i] > 1, ratio[i], -1/ratio[i])
}

sel$ratio = ratio

write.table(sel,file='FC_alleles_PFAM_all.xls',sep='\t',row.names=F,quote=F)

