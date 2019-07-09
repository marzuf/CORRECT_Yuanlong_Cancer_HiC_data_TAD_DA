# adapted from mlpy Python package

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MZ : RANKS ARE EXPECTED TO BE ZERO-BASED !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

harm <- function(n) {
  h  <- 0
  for(i in seq_len(n))
    h <- h+1/i
  return(h)
}



e_harm <- function(n) {
  return(0.5 * harm(floor(n / 2.0)))
}


o_harm <- function(n) {
  return(harm(n) - 0.5 * harm(floor(n / 2.0)))
}

a_harm <- function(n){
  return( ifelse(n %% 2, o_harm(n), e_harm(n)))
}



canberra <- function(x, y){
    return( sum(abs(x-y)/(abs(x) + abs(y))) )
}

    


canberra_location <- function(x, y, k=NULL){
#    """Returns the Canberra distance between two position lists,
#    `x` and `y`. A position list of length P contains the position 
#    (from 0 to P-1) of P elements. k is the location parameter,
#    if k=None will be set to P.
#    """

    if(length(x) != length(y))
       stop("x, y: shape mismatch")
    
    if(is.null(k))
        k = length(x)

    if(k <= 0 | k > length(x)) {
        stop(paste0("k must be in [1, ", length(x), "\n"))
    }
    xx <- ifelse(x+1 < k+1, x+1, k+1)
    yy <- ifelse(y+1 < k+1, y+1, k+1)    


    d=sum(abs(xx-yy)/(xx+yy))
  
    return(d)
}


canberra_stability <- function(x, k=NULL, checkZero=TRUE){
#    """Returns the Canberra stability indicator between N position
#    lists, where `x` is an (N, P) matrix. A position list of length 
#    P contains the position (from 0 to P-1) of P elements. k is 
#    the location parameter, if k=None will be set to P. The lower 
#    the indicator value, the higher the stability of the lists.

#    The stability is computed by the mean distance of all the 
#    (N(N-1))/2 non trivial values of the distance matrix (computed
#    by canberra_location()) scaled by the expected (average) 
#    value of the Canberra metric.

#    Example:
  
  # x = matrix(c(2,4,1,3,0,3,4,1,2,0,2,4,3,0,1), byrow = T, ncol=5)
  # 
  # canberra_stability(x,3)
#    0.74862979571499755
#    """
    if(checkZero) {
      if(! any(x == 0)) {
        stop("! ranking should be 0-based ! \n")
      }
    }
  
    if(is.null(k))
        k = ncol(x)

    if(k <= 0 | k > ncol(x))
        stop(paste0("k must be in [1, ", ncol(x), "]\n"))

    d <- 0.0

    for(i in 1:(nrow(x)-1)){
        for(j in (i+1):nrow(x)){
            d <- d+canberra_location(x[i,], x[j,], k)
        }
    }

    expected = canberra_expected(ncol(x), k)

    return( (d / ((nrow(x)*(nrow(x)-1)) / 2.0) )/expected )
}


canberra_location_expected <- function(p, k=NULL){
#    """Returns the expected value of the Canberra location distance,
#    where `p` is the number of elements and `k` is the number of 
#    positions to consider.
#    """
    if(is.null(k)){
        k <- p
    }
    if(k <= 0 | k > p){
        stop(paste0("k must be in [1, %i]",p, "\n"))
    }
    return(canberra_expected(p,k))
}


canberra_expected <- function(n, k){
    
    mysum <- 0
    
    for(t in 1:k) {
      mysum <- mysum + t * (a_harm(2 * k - t) - a_harm(t));
    }
    # as implemented, the <>harm functions take only 1 integer !!! cannot vectorize !
    # mysum <- sum( tvect * (a_harm(2 * k - tvect) - a_harm(tvect)) )
    # return(2.0 / n * mysum + (2.0 * (n - k) / n) * (2 * (k + 1) * (harm(2 * k + 1) - harm(k + 1)) - k))
    
    
    d = 2.0 / n * mysum + (2.0 * (n - k) / n) * (2 * (k + 1) * (harm(2 * k + 1) - harm(k + 1)) - k)
  
    return(d)
    
    
}

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  RANKS ARE EXPECTED TO BE ZERO-BASED !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



x = matrix(c(2,4,1,3,0,3,4,1,2,0,2,4,3,0,1), byrow = T, ncol=5)

canberra_stability(x,3)
# [1] 0.7486298
canberra_stability(x)
# [1] 0.6688182
canberra_stability(x) == canberra_stability(x,5) 

