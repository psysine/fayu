#!/usr/bin/Rscript
raw1 = read.table('raw', header=TRUE)
attach(raw1)
raw = data.frame(knew=knew, know=know, sincelast=now-pmax(lastknow,lastdunno), a=lastknow-lastdunno-(1-!!lastdunno)*(firstseen-400), b=b, firstseen=firstseen, now=now, lastdunno=lastdunno, lastknow=lastknow, nseen=nseen, nknow=nknow)
detach(raw1)

raw$lsl = log(1+raw$sincelast)
raw$la = log(1+pmax(raw$a,0))
raw$lb = log(1+raw$b)
raw$c = raw$nknow/raw$nseen

formulae=c(know ~ lsl,
           know ~ lsl + I(la*lsl) + I(la^2*lsl) + I(lb*lsl) + I(c*lsl) + la + I(la^2) + lb + c,
           know ~ lsl + c,
           know ~ lsl + I(lb*lsl) + I(c*lsl) + lb + c)

divide = function(d) list(d$knew & !d$b,
                          d$knew & d$b,
                          !d$knew & !d$b,
                          !d$knew & d$b)

sections = c("knew, b=0","knew, b>0","!knew, b=0","!knew, b>0")
div=divide(raw)
for(j in 1:4) {
  rawj=raw[div[[j]],]
  m = glm(formulae[[j]], data=rawj, family='binomial')
  cat(sections[j], "\n")
  cat(names(m$coef), "\n", sep="   ")
  cat(m$coef, "\n", sep="   ")
  }
