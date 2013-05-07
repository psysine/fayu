#!/usr/bin/Rscript
#knew know firstseen lastknow lastdunno nknow ndunno maxkt max2kt sumkt sumdt sumlkt sumldt now k1 k2 k3 k4 k5 k6 k7 k8 k9 k10
raw = read.table('raw', header=TRUE)

raw$sincelast = raw$now-pmax(raw$lastknow, raw$lastdunno)
raw$a=raw$lastknow-ifelse(raw$lastdunno>0, raw$lastdunno, raw$firstseen-400)

raw$lsl = log(1+raw$sincelast)
raw$la = log(1+pmax(raw$a,0))

raw$nseen = raw$nknow+raw$ndunno
raw$c = raw$nknow/raw$nseen
raw$d = raw$nknow/(raw$ndunno+1)
raw$lnseen = log(raw$nseen)
raw$lmaxkt = log(1+raw$maxkt)
raw$lmax2kt = log(1+raw$max2kt)
raw$lsumkt = log(1+raw$sumkt)
raw$lsumdt = log(1+raw$sumdt)
raw$lsumlkt = log(1+raw$sumlkt)
raw$lsumldt = log(1+raw$sumldt)

formulae=c(know ~ lsl,
           know ~ lsl + I(la*lsl)+la + I(la^2*lsl)+I(la^2) + I(lmaxkt*lsl)+lmaxkt + I(lmax2kt*lsl)+lmax2kt + I(d*lsl)+d + I(lsumkt*lsl)+lsumkt + I(lsumdt*lsl)+lsumdt + I(sumlkt*lsl)+sumlkt + I(sumldt*lsl)+sumldt  + I(lsumldt*lsl)+lsumldt + I(c*lsl)+c + I((k2+k3+k4)*lsl)+I(k2+k3+k4) + I((k5+k6+k7)*lsl)+I(k5+k6+k7) + I(lnseen*lsl)+lnseen,
           know ~ lsl + c,
           know ~ lsl + I(lmaxkt*lsl)+lmaxkt + I(lmax2kt*lsl)+lmax2kt + I(d*lsl)+d + I(lsumkt*lsl)+lsumkt + I(lsumdt*lsl)+lsumdt + I(sumlkt*lsl)+sumlkt + I(sumldt*lsl)+sumldt  + I(lsumldt*lsl)+lsumldt + I(c*lsl)+c + I((k2+k3+k4)*lsl)+I(k2+k3+k4) +I((k5+k6+k7)*lsl)+I(k5+k6+k7) + I(lnseen*lsl)+lnseen)

divide = function(d) list(d$knew & !d$maxkt,
                          d$knew & d$maxkt,
                          !d$knew & !d$maxkt,
                          !d$knew & d$maxkt)
if(1) {
  sections = c("knew,maxkt=0","knew,maxkt>0","!knew,maxkt=0","!knew,maxkt>0")
  div=divide(raw)
  for(j in 1:4) {
    rawj=raw[div[[j]],]
    m = glm(formulae[[j]], data=rawj, family='binomial')
    cat(sections[j], "\n")
    cat(length(m$coef), "\n")
    cat(gsub(' ', '', names(m$coef)), "\n", sep="   ")
    co = m$coef
    co[is.na(co)] = 0
    cat(co, "\n", sep="   ")
    }
} else {
  div=divide(raw)
  raw1=raw[div[[1]],]
  m1 = glm(formulae[[1]], data=raw1, family='binomial')
  raw2=raw[div[[2]],]
  m2 = glm(formulae[[2]], data=raw2, family='binomial')
  raw3=raw[div[[3]],]
  m3 = glm(formulae[[3]], data=raw3, family='binomial')
  raw4=raw[div[[4]],]
  m4 = glm(formulae[[4]], data=raw4, family='binomial')
  m=list(m1,m2,m3,m4)
  print(sum(sapply(m, function(x) x$deviance)))
  }
