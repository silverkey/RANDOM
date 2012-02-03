# This script prepare the data.frame with the coordinates of all the repeats.
# Probably sooner or later this could be dane more easily taking advantage of the
# GenomicFeatures package, but at the moment it gives an error probably due to the
# big amount of data to download. So it is better first to download the repeats tables
# and then run this script.

# You have to download UCSC tables containing repeats mapping with the following command:
# wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/*rmsk*'
# Then run this script to build the data.frame to be used in other scripts.

# ----------------- #
#   RUNTIME VARS    #
# ----------------- #
# The directory in which repeats table have been downloaded
rep.dir = '/Users/remo/Desktop/RIKEN/UCSC/REPEATS'
# The name you want to give to the file containing the repeats coordinates
outname = 'rmsk.mm9.RData'

# ----------------- #
#     FUNCTIONS     #
# ----------------- #
get.header = function(tabs) {
  header = readLines(tabs[grep('chr.+sql',tabs)][1])
  header = header[grep('  `',header)]
  header = sub('  `(.+)`.+,',header,replacement="\\1")
  header
}

make.rep.table = function(data.file,header) {
  rmsk = matrix(nrow=0,ncol=7)
  colnames(rmsk) = c('genoName','genoStart','genoEnd','strand','repName','repClass','repFamily')
  for(i in 1:length(data.file)) {
    tab = read.delim(data.file[i],header=F)
    colnames(tab) = header
    rmsk = rbind(rmsk,tab[,colnames(rmsk)])
    message(paste('Loaded table ',data.file[i],'....',sep=''))
  }
  rmsk
}

# ----------------- #
#      SCRIPT       #
# ----------------- #
setwd(rep.dir)
tabs = dir()
data.file = tabs[grep('chr.+txt',tabs)]
header = get.header(tabs)
rmsk = make.rep.table(data.file,header)
save(rmsk,file=outname)
