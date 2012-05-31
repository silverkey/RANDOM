# ----------------- #
#   RUNTIME VARS    #
# ----------------- #
# The directory and the file inside it to annotate
sel.dir = '/Users/remo/Desktop/RIKEN/FINAL_SELECTED_DE/MAPPING_ANALYSIS/SECOND/'
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
calculate.range.overlap = function(ranges,features,r.dir) {
  f.ov2 = findOverlaps(ranges,features,select='first')
  rows = 1:length(f.ov2)
  f.ov = data.frame(row.names=rows)
  f.ov[,1] = rows
  f.ov[,2] = f.ov2
  f.ov = na.omit(f.ov)
  id = as.character(elementMetadata(ranges[f.ov[,1]])$id)
  r.name = as.character(elementMetadata(features[f.ov[,2]])$r.name)
  r.class = as.character(elementMetadata(features[f.ov[,2]])$r.class)
  r.family = as.character(elementMetadata(features[f.ov[,2]])$r.family)
  r.strand = as.character(strand(features[f.ov[,2]]))
  res = as.data.frame(unique(cbind(id,r.name,r.class,r.family,r.strand,r.dir)))
  no.res = ranges
  if(ncol(res)>1) no.res = ranges[-f.ov[,1]]
  list(res=res,no.res=no.res)
}

antisense.strand = function(strand) {
  a = 'NA'
  if(strand == '+') a = '-'
  if(strand == '-') a = '+'
  if(strand == '1') a = '-1'
  if(strand == '-1') a = '1'
  if(a == 'NA') stop(paste('Cannot recognize strand: ',strand,sep=''))
  a
}

granges.antisense = function(ranges) {
  strand = as.data.frame(ranges)$strand
  strand(ranges)=Rle(unlist(lapply(strand,antisense.strand)))
  ranges
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
ncol.f = ncol(source.tab)
selected = read.delim(file=data)[,c(id,chr,start,end,strand)]
colnames(selected) = c('id','chr','start','end','strand')
#selected$strand = '*';

features.res.cn = c('id','r.name','r.class','r.family','r.strand','r.dir')
features.res = matrix(ncol=length(features.res.cn),nrow=0)
features.tab = matrix(ncol=ncol.f+length(features.res.cn)-1,nrow=0)
colnames(features.res) = features.res.cn

# Build the ranges of the fragments you want to test the overlap with features
ranges = GRanges(seqnames=Rle(selected$chr),
                 ranges=IRanges(start=selected$start,end=selected$end),
                 strand=Rle(selected$strand),
                 id=selected$id)
  
lres = calculate.range.overlap(ranges,repeats,'sense')
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

ranges = granges.antisense(ranges)

lres = calculate.range.overlap(ranges,repeats,'antisense')
ranges = lres$no.res
if(ncol(lres$res)==ncol(features.res)) features.res = rbind(features.res,lres$res)

no.rep = as.data.frame(lres$no.res)
if(nrow(no.rep)>=1) {
  no.rep$r.name=NA
  no.rep$r.class=NA
  no.rep$r.family=NA
  no.rep$r.strand=NA
  no.rep$r.dir=NA
  features.res = rbind(features.res,no.rep[,c('id','r.name','r.class','r.family','r.strand','r.dir')])
}

features.tab = rbind(features.tab,merge(source.tab,features.res,by.x=id,by.y='id'))

write.table(unique(features.tab),file='RESULTS_FEATURE_REPEAT_ANNOTATED.xls',sep="\t",row.names=F,quote=F)

# CHECK IF THE ANTISENSE ARE ENRICHED IN A SPECIFIC LOCATION
sense=table(features.tab[features.tab$r.dir=='sense',]$r.class)
anti=table(features.tab[features.tab$r.dir=='antisense',]$r.class)
tot.sense=sum(table(features.tab[features.tab$r.dir=='sense',]$r.class))
tot.anti=sum(table(features.tab[features.tab$r.dir=='antisense',]$r.class))
c=cbind(sense,anti,tot.sense,tot.anti)
test = function(x) { prop.test(c(x[1],x[2]),c(x[3],x[4]),alternative='l')$p.value }
apply(c,1,test)
par(mfrow=c(1,2))
barplot(sense[sense>=20]/tot.sense*100,ylim=c(0,40),main=paste('Percentage DE clusters overlapping repeats in sense orientation (n=',tot.sense,')',sep=''))
barplot(anti[anti>=20]/tot.anti*100,ylim=c(0,40),main=paste('Percentage DE clusters overlapping repeats in antisense orientation (n=',tot.anti,')',sep=''))
