distmatrix=function(coords){

  if(sum(is.na(coords)) > 0) stop("There are some NA values in the specified coordinates")
  if (!is.numeric(coords)& !is.data.frame(coords)) stop(" 2D coordinates must be a numeric matrix or data.frame")
  if (!is.matrix(coords)) coords=as.matrix(coords)
  if (ncol(as.matrix(coords)) !=2) stop("2D coordinates must be specified")

  out=dist(coords)

  return(list(distmatrix=out$dist))

}