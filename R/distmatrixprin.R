####4.Distances matrix############################

dist<-function(coords){

  n=nrow(coords)
  dist=matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      dist[i,j] = sqrt(sum((coords[i,]-coords[j,])^2));
    }
  }
  dist=dist+t(dist)/2;
  return(list(dist=dist))
}
