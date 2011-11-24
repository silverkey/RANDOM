cal.cv = function(m) {
  cs = c()
  r = 3
  v = c(1:r)
  n = length(m)
  for(i in 1:(n/r)){
    s = as.numeric(m[v])
    c = sd(s)/mean(s)
    v = v + r
    cs = c(cs,c)
  }
  cs
}

filt.cv = function(row) {
  cutoff = 0.8
  res = c()
  vec = na.omit(row)
  r = sum(vec>cutoff)
  res = c(res,r)
  res
}

CVfile = 'cv.pdf'
CVtitle = 'cv'

raw = read.delim(file='raw.txt')
cv.m = apply(raw[,-1],1,cal.cv)
cv.m = t(cv.m)
raw.cv = as.data.frame(cv.m)
cv.calc = apply(raw.cv,1,filt.cv)
cv.selected = as.vector(raw[cv.calc==0,1])

pdf(file=CVfile,paper='a4r',width=8.3,height=11.7,pointsize=8,title=CVtitle)
vec = as.vector(unlist(na.omit(raw.cv)))
hist(vec,col=8,breaks=seq(0,2,0.01),xlab='coefficient of variation',ylab='number of clusters',main='Coefficients of Variation Raw Level 2 Data')
dev.off()

