#!/usr/bin/R

# Script by Remo Sanges 21/05/2012

# Finally GOstat on neglected organisms!!!

# Write the table in output
outname = 'up_GO_significant.xls'

# The file of the selected genes. Must contains a single column with
# the probe id of the differentially expressed genes
selected.filename = 'up.txt'

# The file containing the association of GO id to the probe id.
# It must contains 2 columns, the first with the probe id and 
# the second with the associated GO id.
go.annotation.filename = 'GO_array_annotation.txt'

# The file with the definition of the GO id. It contains 2 columns
# the first containing the GO id and the second with the definitions
# This file come out using the script parse_standard_go_table.pl
# on the official GO file GO.terms_alt_ids
go.definition.filename = 'go_definition.txt'

# The cutoff on the adjusted p.value
pcutoff = 0.1

# The cutoff on the minimum number of probes selected belonging to the class
ncutoff = 10

# The function to test GO enrichments in every GO division
testGoDivision = function(ass,sel,def,classname,ncutoff,pcutoff) {

  res = NULL
  
  ass = ass[ass[,2] %in% def[,1],]
  sel.ass = ass[ass[,1] %in% sel[,1],]
  classes = as.vector(unique(sel.ass[,2]))
  tot.uni = length(as.vector(unique(ass[,1])))
  tot.sel = length(as.vector(unique(sel.ass[,1])))

  for(i in 1:length(classes)) {
    
    class = classes[i]
    
    n.uni = length(as.vector(unique(ass[ass[,2]==class,1])))
    n.sel = length(as.vector(unique(sel.ass[sel.ass[,2]==class,1])))
    
    p.uni = n.uni/tot.uni*100
    p.sel = n.sel/tot.sel*100
    
    if(n.sel >= ncutoff) {
      if(p.sel >= p.uni) {
        probes = paste(sel.ass[sel.ass[,2]==class,1],collapse=',')
        class.def = as.character(def[def[,1]==class,2])
        p.value = prop.test(c(n.sel,n.uni),c(tot.sel,tot.uni),alternative='g')$p.value
        res = rbind(res,data.frame(GO.id=class,n.uni=n.uni,n.sel=n.sel,p.uni=p.uni,p.sel=p.sel,p.value=p.value,adj.p=1,definition=class.def,probes=probes))
        res = na.omit(res)
      }
    }
  }
  adj.pval = p.adjust(res$p.value)
  res$adj.p = adj.pval
  
  sig = res[res$adj.p<=pcutoff,]
  sig = sig[order(sig$p.value),]
  
  write.table(sig,file=paste(classname,ncutoff,outname,sep='_'),sep="\t",quote=F,row.names=F)
}

sel = read.table(file=selected.filename,sep="\t",header=F,quote="")
ass = read.table(file=go.annotation.filename,sep="\t",header=F,quote="")
def = read.table(file=go.definition.filename,sep="\t",header=F,quote="")
colnames(def) = c('goid','goname','class')
def.bp = def[def$class == 'P',]
def.mf = def[def$class == 'F',]
def.cc = def[def$class == 'C',]

testGoDivision(ass,sel,def.bp,'BP',ncutoff,pcutoff)
testGoDivision(ass,sel,def.mf,'MF',ncutoff,pcutoff)
testGoDivision(ass,sel,def.cc,'CC',ncutoff,pcutoff)
