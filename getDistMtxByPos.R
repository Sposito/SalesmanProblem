#' Returns a a distance matrix.
#' 
#' Given a matrix where the rows are points and the columns are coordinates
#' it returns another matrix that computes the distances between the points
getDistMtxByPos <- function(pos){
    rows = nrow(pos)
    distances = matrix( rep(0,rows ** 2), nrow=rows, ncol=rows, byrow = FALSE)
  
    for (i in 1:rows){
        for (j in 1:rows){
            xI <- pos[i, 1];
            yI <- pos[i, 2];
            xJ <- pos[j, 1];
            yJ <- pos[j, 2];
      
            distances[i,j] = sqrt( (xJ - xI) ** 2 + (yJ - yI) ** 2);
        } 
    }
    return (distances);
}

# myPos = matrix( c(33,  1,  0,  0,  2, 30, 31, 29, 31, 32, 30, 31,
#                   30, 31, 29, 31, 32,  0,  1,  0,  0,  2, 30, 31),
#                  nrow=12, ncol=2);
#  
#  dst = getDistMtxByPos(myPos)