# ----------------- #
#   RUNTIME VARS    #
# ----------------- #
# The directory and the file inside it to annotate
sel.dir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/send/sel_coord'
data = 'RESULTS_FEATURE_ANNOTATED.xls'
# The name with path of the repeat file created with the script build_repeats_matrix.R
dbdir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/DATA'
repfile = paste(dbdir,'rmsk.mm9.RData',sep='/')
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
calculate.range.overlap = function(ranges,features) {
  f.ov = as.matrix(findOverlaps(ranges,features))
  id = as.character(elementMetadata(ranges[f.ov[,1]])$id)
  r.name = as.character(elementMetadata(features[f.ov[,2]])$r.name)
  r.class = as.character(elementMetadata(features[f.ov[,2]])$r.class)
  r.family = as.character(elementMetadata(features[f.ov[,2]])$r.family)
  r.strand = as.character(strand(features[f.ov[,2]]))
  res = as.data.frame(unique(cbind(id,r.name,r.class,r.family,r.strand)))
  no.res = ranges
  if(ncol(res)>1) no.res = ranges[-f.ov[,1]]
  list(res=res,no.res=no.res)
}

# ----------------- #
#      SCRIPT       #
# ----------------- #
library("GenomicFeatures")
setwd(sel.dir)

load(repfile)
repeats = GRanges(seqnames=Rle(rmsk$genoName),
                  ranges=IRanges(start=rmsk$genoStart,end=rmsk$genoEnd),
                  strand=Rle(rmsk$strand),
                  r.name=rmsk$repName,r.class=rmsk$repClass,r.family=rmsk$repFamily)

source.tab = read.delim(file=data)
selected = read.delim(file=data)[,c(id,chr,start,end,strand)]
colnames(selected) = c('id','chr','start','end','strand')
selected$strand = '*';

ncol.f = ncol(source.tab)
features.res = matrix(ncol=5,nrow=0)
features.tab = matrix(ncol=ncol.f+5-1,nrow=0)
colnames(features.res) = c('id','r.name','r.class','r.family','r.strand')

# Build the ranges of the fragments you want to test the overlap with features
ranges = GRanges(seqnames=Rle(selected$chr),
                 ranges=IRanges(start=selected$start,end=selected$end),
                 strand=Rle(selected$strand),
                 id=selected$id)
  
lres = calculate.range.overlap(ranges,repeats)
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

no.rep = as.data.frame(lres$no.res)
if(nrow(no.rep)>=1) {
  no.rep$r.name=NA
  no.rep$r.class=NA
  no.rep$r.family=NA
  no.rep$r.strand=NA
  features.res = rbind(features.res,no.rep[,c('id','r.name','r.class','r.family','r.strand')])
}

features.tab = rbind(features.tab,merge(source.tab,features.res,by.x=id,by.y='id'))

write.table(unique(features.tab),file='RESULTS_FEATURE_REPEAT_ANNOTATED.xls',sep="\t",row.names=F,quote=F)
