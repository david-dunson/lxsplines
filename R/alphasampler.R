#######################################################################
# function sampA
#   input: Y - data
#          w - X matrix
#          tau - estimated covariance
#          E   - prior expectation
#          V   - prior variance
#          amin - maximum value of a
#          amax - minimum value of a
#          omega - Weight matrix
#    ouptput: Gibbs sample of a
#########################################################################
sampA<-function(Y,w,tau,E,V,amin,amax,omega){
#Sample the A matrix given a Y and a precomputed W matrix

  tW = w*matrix(omega,nrow=nrow(w),ncol=ncol(w))
  EV = 1/(tau*(t(tW)%*%w + 1/V)); #posterior mean   
  ES = EV*(tau*t(tW)%*%Y + E/V);  #posterior variance
  UB = pnorm(1,ES,sqrt(EV)); 
  LB = pnorm(-1,ES,sqrt(EV)); 
  alpha = qnorm(runif(1,LB,UB),ES,sqrt(EV));           
  return(alpha);                                             
}

