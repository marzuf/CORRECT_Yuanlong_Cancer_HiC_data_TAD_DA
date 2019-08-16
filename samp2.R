
require(foreach)


myfunc <- function(n) {
    set.seed(123)
    myvector <- c("a", "b", "c", "d")
    all_x <- foreach(i=1:n, .combine='c') %do% {       
        mysubfunc(myvector)           
    }
    all_x
}



mysubfunc <- function(myvector) {
    sample(myvector, size=length(myvector), replace=F)     
}




