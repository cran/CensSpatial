atypicalprin=function(w){
  res=w$respper
  a=res[res[,1]=="atypical obs",]

  sm=w$smper
  b=sm[sm[,1]=="atypical obs",]

  ev=w$expvper
  c=ev[ev[,1]=="atypical obs",]
  s=list(RP=a,SP=b,EP=c)
  return(s)
}
