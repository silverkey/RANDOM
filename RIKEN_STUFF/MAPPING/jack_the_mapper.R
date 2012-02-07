# This script build the genomic annotations by using the GenomicFeatures package
# It is very fast, but you cannot save yet the structures. Therefore after saving the
# sqlite db you will have to reload it and to recalculate each time all the features
# (transcripts, promoters, exons, utr etc.). However this is not taking very long time....

# ----------------- #
#   JUST BEGINS..   #
# ----------------- #
# Lazy loading of a map table between ensembl ids and symbols
# because these info do not seem to be loaded by GenomicFeatures
# This will create into the workspace a data.frame variable named
# "map" with the pairings ensembl.gene.id <------> symbol
source('make.id.map.R')

# ----------------- #
#   RUNTIME VARS    #
# ----------------- #
# The directory and the file inside it to annotate
sel.dir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/send/sel_coord'
data = 'RESULTS_COORD.xls'
# The directory containing the database or in which you will create it
dbdir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/DATA'
genome = 'mm9'
tablename = 'ensGene'
download = 'F'
# Length of promoters
proximal = 1000
# Columns containing the postional info in the files to annotate
# All the next fields are mandatory and id has to be unique!
id = 1
chr = 7
start = 8
end = 9
strand = 10

# ----------------- #
#     FUNCTIONS     #
# ----------------- #
get.promoters = function(transcripts,proximal) {
  pos = strand(transcripts) == '+'
  promoter.start = as.vector(ifelse(pos,start(transcripts)-proximal,end(transcripts)+1))
  promoter.end = as.vector(ifelse(pos,start(transcripts)-1,end(transcripts)+proximal))
  promoter.strand = ifelse(pos,'+','-')
  promoter.chr = seqnames(transcripts)
  promoter.tx_name = unlist(elementMetadata(transcripts[,'tx_name'])[,1])
  promoter.gene_id = unlist(elementMetadata(transcripts[,'gene_id'])[,1])
  promoters = GRanges(seqnames=promoter.chr,
                      ranges=IRanges(start=promoter.start,end=promoter.end),
                      strand=promoter.strand,
                      tx_name=promoter.tx_name,gene_id=promoter.gene_id)
  promoters
}

calculate.range.overlap = function(ranges,features,f.name) {
  f.ov = as.matrix(findOverlaps(ranges,features))
  id = as.character(elementMetadata(ranges[f.ov[,1]])$id)
  f.gene = as.character(elementMetadata(features[f.ov[,2]])$gene_id)
  f.strand = as.character(strand(features[f.ov[,2]]))
  res = as.data.frame(unique(cbind(id,f.name,f.gene,f.strand)))
  if(ncol(res)<4) no.res = ranges
  if(ncol(res)==4) no.res = ranges[-f.ov[,1]]
  list(res=res,no.res=no.res)
}

calculate.rangelist.overlap = function(ranges,features,f.name,txid.gid) {
  f.ov = as.matrix(findOverlaps(ranges,features))
  id = as.character(elementMetadata(ranges[f.ov[,1]])$id)
  f.gene = as.character(txid.gid[names(features[f.ov[,2]]),2])
  f.strand = as.character(txid.gid[names(features[f.ov[,2]]),3])
  res = as.data.frame(unique(cbind(id,f.name,f.gene,f.strand)))
  if(ncol(res)<4) no.res = ranges
  if(ncol(res)==4) no.res = ranges[-f.ov[,1]]
  list(res=res,no.res=no.res)
}

# ----------------- #
#      SCRIPT       #
# ----------------- #
library("GenomicFeatures")

if(download == 'T') {
  transdb = makeTranscriptDbFromUCSC(genome=genome,tablename=tablename)
  saveFeatures(transdb,file=paste(dbdir,'/',genome,'.',tablename,'.','sqlite',sep=''))
} else {
  transdb = loadFeatures(file=paste(dbdir,'/',genome,'.',tablename,'.','sqlite',sep=''))
}

setwd(sel.dir)

# Completely noncodings transcripts will overlap exons but not cds nor utr
# We use a hierarchy order for which every time a range is overlapping a feature than
# it is associated to the feature and cut-out from the ranges, so that it cannot overlapping
# with the next order feature.
# Order to test overlap:
# 1) promoters
# 2) utr
# 3) cds
# 4) exon
# 5) intron

# Build Features
transcripts = transcripts(transdb,columns=c("tx_id","tx_name","gene_id")) # range
cds = cdsBy(transdb,by='tx',use.names=T)                                  # rangelist
exons = exonsBy(transdb,by='tx',use.names=T)                              # rangelist
introns = intronsByTranscript(transdb,use.names=T)                        # rangelist
utr5 = fiveUTRsByTranscript(transdb,use.names=T)                          # rangelist
utr3 = threeUTRsByTranscript(transdb,use.names=T)                         # rangelist
promoters = get.promoters(transcripts,proximal)                           # range

# Build a comfortable table to associate transcripts->genes->strands
txid.gid = as.data.frame(cbind(
                         unlist(elementMetadata(transcripts[,'tx_name'])[,1]),
                         unlist(elementMetadata(transcripts[,'gene_id'])[,1]),
                         as.character(strand(transcripts))))
colnames(txid.gid) = c('tx_id','gene_id','strand')
rownames(txid.gid)=txid.gid$tx_id

source.tab = read.delim(file=data)
ncol.f = ncol(source.tab)
selected = read.delim(file=data)[,c(id,chr,start,end,strand)]
colnames(selected) = c('id','chr','start','end','strand')
selected$strand = '*';

features.res.cn = c('id','f.name','f.gene','f.strand')
features.res = matrix(ncol=length(features.res.cn),nrow=0)
features.tab = matrix(ncol=ncol.f+length(features.res.cn)-1,nrow=0)
colnames(features.res) = features.res.cn

# Build the ranges of the fragments you want to test the overlap with features
ranges = GRanges(seqnames=Rle(selected$chr),
                 ranges=IRanges(start=selected$start,end=selected$end),
                 strand=Rle(selected$strand),
                 id=selected$id)
  
lres = calculate.range.overlap(ranges,promoters,'promoter')
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

lres = calculate.rangelist.overlap(ranges,utr5,'utr5',txid.gid)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

lres = calculate.rangelist.overlap(ranges,utr3,'utr3',txid.gid)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

lres = calculate.rangelist.overlap(ranges,cds,'cds',txid.gid)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

lres = calculate.rangelist.overlap(ranges,exons,'exon',txid.gid)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

lres = calculate.rangelist.overlap(ranges,introns,'intron',txid.gid)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

intergenic = as.data.frame(lres$no.res)
if(nrow(intergenic)>=1) {
  intergenic$f.name='intergenic'
  intergenic$f.gene=NA
  intergenic$f.strand=NA
  features.res = rbind(features.res,intergenic[,features.res.cn])
}

features.tab = rbind(features.tab,merge(source.tab,features.res,by.x=id,by.y='id',all.x=T))
col.order = c(colnames(features.tab),'sym')
features.tab = merge(features.tab,map,by.x='f.gene',by.y='ens',all.x=T)
features.tab = features.tab[,col.order]
write.table(unique(features.tab),file='RESULTS_FEATURE_ANNOTATED.xls',sep="\t",row.names=F,quote=F)
