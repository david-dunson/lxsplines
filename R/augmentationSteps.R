################################################################
# augmentFunctions: 
# Purpose: Creates the polya Gamma Augment
#
#
#
#
#
################################################################


#Samples from a polya gamma PG(a,b) distribution
.neg_binomial_augment<-function(y,r_disp,vecMean,truncation = 350){
        temp.v <- matrix(c(y+r_disp,vecMean),nrow=2,byrow=TRUE)
        
        
        omega  <-  polyaGamma(temp.v , truncation)
    
        return(t(omega))
}

#samples iterativly from a poly gamma PG(1,b) distribution
#for more efficient samples
.logig_binomial_augment<-function(y,vecMean){
  #as a temporary measure do the random sum 
  if (ncol(y) != 2) #not enough enformation
  { stop("Binomial not properly described!")}
  temp.v <- t(y)
  temp.v[1,] <- temp.v[2,]
  temp.v[2,] <- vecMean
  #truncate temporarily
  omega  <-  polyaGamma(temp.v , 1000)
  return(t(omega)); 
};