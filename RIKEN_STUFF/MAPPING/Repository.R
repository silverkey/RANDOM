# ----------------- #
#   RUNTIME VARS    #
# ----------------- #

# ----------------- #
#     FUNCTIONS     #
# ----------------- #

# ----------------- #
#      SCRIPT       #
# ----------------- #

# ----------------- #
#    REPOSITORY     #
# ----------------- #
antisense.strand = function(strand) {
  a = 'NA'
  if(strand == '+') a = '-'
  if(strand == '-') a = '+'
  if(strand == '1') a = '-1'
  if(strand == '-1') a = '1'
  if(a == 'NA') stop(paste('Cannot recognize strand: ',strand,sep=''))
  a
}
