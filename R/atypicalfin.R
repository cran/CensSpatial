atypical=function(w){
  if (class(w)!="localinfmeas") stop("an object of the class localinfmes must be provided")
 out=atypicalprin(w)
   return(out)
}
