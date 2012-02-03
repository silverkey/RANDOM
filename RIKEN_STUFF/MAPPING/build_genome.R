# ----------------- #
#   RUNTIME VARS    #
# ----------------- #
genome = 'mm9'
tablename = 'ensGene'
download = 'F'
proximal = 1000

# ----------------- #
#     FUNCTIONS     #
# ----------------- #
get.promoters = function(transcripts,proximal) {
  pos = strand(transcripts) == '+'
  promoter.start = as.vector(ifelse(pos,start(transcripts)-proximal,end(transcripts)))
  promoter.end = as.vector(ifelse(pos,start(transcripts),end(transcripts)+proximal))
  promoter.strand = ifelse(pos,'+','-')
  promoter.chr = seqnames(transcripts)
  #promoter.id = elementMetadata(transcripts[,'tx_name'])[,1]
  promoters = GRanges(seqnames=promoter.chr,ranges=IRanges(start=promoter.start,end=promoter.end),strand=promoter.strand)
  promoters
}

# ----------------- #
#      SCRIPT       #
# ----------------- #
library("GenomicFeatures")

if(download == 'T') {
  transdb = makeTranscriptDbFromUCSC(genome=genome,tablename=tablename)
  saveFeatures(transdb, file=paste(genome,tablename,'sqlite',sep='.'))
} else {
  transdb = loadFeatures(file=paste(genome,tablename,'sqlite',sep='.'))
}

# Completely noncodings transcripts will overlap exons but not cds nor utr
# Order to test overlap:
# -promoters
# -utr
# -cds
# -exon
# -check for noncodingness
# -intron
# -check for intergenicness

transcripts = transcripts(transdb)
cds = cdsBy(transdb,by='tx')
exons = exonsBy(transdb,by='tx')
introns = intronsByTranscript(transdb)
utr5 = fiveUTRsByTranscript(transdb)
utr3 = threeUTRsByTranscript(transdb)
promoters = get.promoters(transcripts,proximal)
