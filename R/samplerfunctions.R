###############################################################################
#
#
#
#
#
#
###############################################################################
lxfit  <- function(x,y, anal_type = 'normal', isIncrease = TRUE,mEXTREMA = 2,
                   nsamps = 50000,
                   nburn =  10000,
                   a_list = c(1,1,1,1.15,1.3,1.7,2,2.6,2.0,3.5,5,9,12,24,30,50)){
  #make sure the y variable is a matrix
  t_y = as.matrix(y,ncol=ncol(y))
  #center the x variable to be between 0 and 1
  t_x = matrix((x-min(x))/(max(x)-min(x)),ncol=1 )
  if (nrow(t_x) != nrow(t_y)){
    stop("X and Y vectors must have the same number of rows.")
  }
  
  if (anal_type       == 'normal'){
      analysis = 1
  }else if (anal_type == 'negbin'){
      analysis = 2
  }else if (anal_type == 'binomial'){
      analysis = 3
  }else {
      stop(sprintf("Only Normal 'normal,' Negative Binomial 'negbin', and binomial 'binomial'
                    models are supported at this time. User supplied  %s",analy_type))
  }
  
  
  #error check the JK quantity making sure 1 is in the list and
  #it is sorted
  a_list = sort(a_list)
  if (a_list[1] != 1){
    stop("The minimum value in the anealling list must be 1")
  }
  
  if (length(a_list)==1){
    stop("Anealling list must have more than one element")
  }
  
  if (nsamps < nburn){
    stop("Number of samples must be greater than the number of burnin samples.")
  }
  
  if (nsamps < 0){
    stop("Number of samples must be greater than zero.")
  }
  
  if (nburn < 0){
    warning("Burn in samples less than zero. Setting burn in samples to zero.")
    nburn = 1; 
  }
  
  if (mEXTREMA == 2){
      if (isIncrease == TRUE){
        M = 100
      }else{
        M= -100
      }
    AXTREMA = 2
  }
  if (mEXTREMA == 1)
  {
    if (isIncrease == TRUE){
      M = -100
    }else{
      M=   100
    }
    AXTREMA = 1
  }
  return (.LXsample_linear_2(t_x,t_y,a_list,nsamps,nburn,analysis,M,AXTREMA))
  
}


##############################################################
# Calculate Posterior model probabilities for two changepoints
# nburn - is the number of burnin samples to ignore
# nsamps -  is the maximum sample - needs to be less than 
# the number of samples used in the fit
##############################################################
calculateProbs_2 <- function(fit,nburn,nsamps){
  sa1 = fit$sa1
  sa2 = fit$sa2
  if (nsamps > length(sa2)){
    stop("Number of samples (nsamps) must be less than or equal to the total number of 
          samples")
  }
  SHAPES =  rep(0,length(nburn:nsamps))
  for (ii in 1:length(SHAPES)) {
     tt <- ii + nburn - 1
    #decide 'base shape'
    shape = 1; #'n' shape
    shape =  ((sa1[tt] > -.5)*(sa2[tt]> -.5)*(sa1[tt] < 0.5)*(sa2[tt]< 0.5))*4 + shape; #s shaped
    shape = ((sa1[tt] <= -0.5)*(sa2[tt] >= 0.5)+(sa1[tt] >= 0.5)*(sa2[tt] <= -0.5))*2 + shape# Monotone Decreasing
    shape = ((sa1[tt] <= -0.5)*(sa2[tt] <= -0.5)+(sa1[tt] >= 0.5)*(sa2[tt] >= 0.5))*3 + shape # Monotone Increasing
    shape = ((sa1[tt] > -.5)*(sa1[tt]< 0.5)*(sa2[tt] <= -0.5)+(sa2[tt] > -.5)*(sa2[tt]< 0.5)*(sa1[tt] <= -0.5))*1 +shape#J-shaped
    

    SHAPES[ii] = shape  
  }
  
    pr1 = mean(SHAPES == 5) #S SHAPED
    pr2 = mean(SHAPES == 3) #Monotone Decreasing
    pr3 = mean(SHAPES == 4) #Monotone Increasing
    pr4 = mean(SHAPES == 2) #J shaped
    pr5 = 1-pr1-pr2-pr3-pr4 # n- shaped
  return (list(pr1=pr1,pr2=pr2,pr3=pr3,pr4=pr4,pr5=pr5,SHAPES = SHAPES))
}


#####################################################

#####################################################
##############################################################
# Calculate Posterior model probabilities for two changepoints
# nburn - is the number of burnin samples to ignore
# nsamps -  is the maximum sample - needs to be less than 
# the number of samples used in the fit
##############################################################
calculateProbs_1 <- function(fit,nburn,nsamps){
  sa1 = fit$sa1[nburn:nsamps]
  
  
  pr1 = mean(sa1 < -0.5) #Monotone Increasing
  pr2 = mean(sa1 > 0.5) #Monotone Decreasing
  pr3 = 1-pr1-pr2 #n shaped
  return (list(pr1=pr1,pr2=pr2,pr3=pr3))
}


#####################################################


#####################################################
bayesFactor_2 <- function(HA,HB, out){
  
  if (out$EXTREMA == 2){
      priorprobs <- c(0.3146,
                      0.0963,
                      0.0963,
                      0.2464,
                      0.2464)
      postProbs <-   calculateProbs_2(out,out$mcmcP[1],out$mcmcP[2])
      pprobs <- c(postProbs$pr1,postProbs$pr2,postProbs$pr3,postProbs$pr4,postProbs$pr5)
      
      if (out$isIncreasing){
        testNames <- c("'~'-shaped",
                       "Monotone Decreasing",
                       "Monotone Increasing",
                       "J-shaped",
                       "n-shaped")
          
      }else{
        testNames <- c("Inverse '~'-shaped",
                       "Monotone Increasing",
                       "Monotone Decreasing",
                       "n-shaped",
                       "J-shaped")
      }
      
      if (length(HA) != 5){
        stop("Hyptohesis must have 5 options")
      }
      if (length(HB) != 5){
        stop("Hyptohesis must have 5 options")
      }
      if (sum(HA) >= length(HA)){
        stop("Hypothesis HA must have fewer groups than the total number of available hyptoheses.")
      }
      
      if (sum(HB) >= length(HB)){
        stop("Hypothesis HB must have fewer groups than the total number of available hyptoheses.")
      }
      

  }else{ #When there is just one extrema
    
    priorprobs <- c(0.22,
                    0.22,
                    0.56)
    
    postProbs <-   calculateProbs_1(out,out$mcmcP[1],out$mcmcP[2])
    pprobs <- c(postProbs$pr1,postProbs$pr2,postProbs$pr3)
    
    if (length(HA) != 3){
      stop("Hyptohesis must have 3 options")
    }
    if (length(HB) != 3){
      stop("Hyptohesis must have 3 options")
    }
    
    if (out$isIncreasing){
      testNames <- c("Monotone Decreasing",
                     "Monotone Increasing",
                      "J-shaped")
      
    }else{
      testNames <- c("Monotone Increasing",
                     "Monotone Decreasing",
                     "n-shaped")
    }
    
    
  }
  ###############################################################
  ###############################################################
  cat('\nBayes Factor for\ncomparing Hypotheses:\n\tHA:\n')
  for (ii in 1:length(HA)){
    
    if (HA[ii] > 0){
      cat(sprintf("\t\t%s\n",testNames[ii]))
    }
  }
  
  cat("\t\t-VS-\n\tHB:\n")
  for (ii in 1:length(HB)){
    if (HB[ii] > 0){
      cat(sprintf("\t\t%s\n",testNames[ii]))
    }
  }
  cat("------------------------------------------------------------------------------\n")
  temp = (sum(HB*priorprobs)/sum(HA*priorprobs))
  tempb = (sum(HA*pprobs)/sum(HB*pprobs))
  cat(sprintf("Bayes Factor: %1.3f \n",as.numeric(temp*tempb))) #*(sum(HB*priorprobs)/sum(HA*priorprobs)))
  cat("------------------------------------------------------------------------------\n\n\n") 
  
  cat("------------------------------------------------------------------------------\n")
  cat("Bayes Factor\tStrength of evidence\n")
  cat("-------------------------------------------------------------------------------\n")
  cat("1    to 3.2   	not worth more than a bare mention\n")
  cat("3.2  to 10      positive\n")
  cat("10   to 31.6    strong\n")
  cat("31.6 to 100     very strong\n")
  cat(">100 	        decisive\n")
  cat("\tBased upon H. Jeffreys (1961).")
  
}
  

##############################################################
##############################################################
# Calculate Posterior model probabilities for two changepoints
# nburn - is the number of burnin samples to ignore
# nsamps -  is the maximum sample - needs to be less than 
# the number of samples used in the fit
##############################################################
hcalculateProbs_2 <- function(fit,nburn,nsamps){
  sa1 = fit$sa1
  sa2 = fit$sa2
  if (nsamps > length(sa2)){
    stop("Number of samples (nsamps) must be less than or equal to the total number of 
         samples")
    
  }
  SHAPES =  rep(0,nburn:nsamps)
  for (ii in nburn:nsamps) {
    a.s <- c(sa1[ii],sa2[ii])
    betas = fit$beta_sample[[ii]]
    knots = c(fit$model_sample[[ii]][1,])
    betas = betas[2:length(betas)]
    flat  = rep(0,length(knots)-1)
    n = length(betas)
    #decide 'base shape'
    shape = 1; #'n' shape
    shape =  ((sa1[ii] > -.5)*(sa2[ii]> -.5)*(sa1[ii] < 0.5)*(sa2[ii]< 0.5))*4 + shape; #s shaped
    shape = ((sa1[ii] <= -0.5)*(sa2[ii] >= 0.5)+(sa1[ii] >= 0.5)*(sa2[ii] <= -0.5))*2 + shape# Monotone Decreasing
    shape = ((sa1[ii] <= -0.5)*(sa2[ii] <= -0.5)+(sa1[ii] >= 0.5)*(sa2[ii] >= 0.5))*3 + shape # Monotone Increasing
    shape = ((sa1[ii] > -.5)*(sa1[ii]< 0.5)*(sa2[ii] <= -0.5)+(sa2[ii] > -.5)*(sa2[ii]< 0.5)*(sa1[ii] <= -0.5))*1 +shape#J-shaped
    
    flat[1] = (betas[1] == 0) + (betas[2] ==  0)
    flat[2] = (betas[2] == 0)
    flat[length(flat)] = (betas[n] == 0) + (betas[n-1] == 0)
    flat[length(flat)-1] = (betas[n-1] == 0)
    if (length(flat) > 2){
      for (jj in 3:(n-2)){
        flat[jj-2] = flat[jj-2] + (betas[jj] == 0)
        flat[jj-1] = flat[jj-1] + (betas[jj] == 0)
        flat[jj] = flat[jj] + (betas[jj] == 0)
      }
    }
    intLoc1 = sum(sa1[ii] > knots-0.5)
    intLoc2 = sum(sa2[ii] > knots-0.5)
    temp = sort(c(intLoc1,intLoc2)+1)
    temp = c(1,temp,length(flat)+2)
    flat = c(1,( flat ==  3),1) #if there are three flat betas the region is flat
    
    #different 'checks' depending on the starting shape of the curve  
    if (shape == 5){
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      center = (prod(flat[temp[2]:temp[3]]) ==1)
      right = (prod(flat[temp[3]:temp[4]]) ==1)
      if  (center){
        shape = 4
      }else if (left){
        if  (right){shape = 2}else{shape=3}
      }else  if(right){
        shape = 1    
      }
    }else if(shape == 1 ){
      temp  = unique(temp)
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      right = (prod(flat[temp[2]:temp[3]]) ==1)
      if (left){
        shape  = 3
      }else  if(right){
        shape  = 4    
      }
    }else if(shape == 2){
      temp  = unique(temp)
      left = (prod(flat[temp[1]:temp[2]]) ==1)
      right = (prod(flat[temp[2]:temp[3]]) ==1)
      if (left){
        shape=3
      }else  if(right){
        shape = 2    
      }
    }
    SHAPES[ii] = shape  
  }
  
  pr1 = mean(SHAPES == 5) #S SHAPED
  pr2 = mean(SHAPES == 3) #Monotone Decreasing
  pr3 = mean(SHAPES == 4) #Monotone Increasing
  pr4 = mean(SHAPES == 2) #J shaped
  pr5 = 1-pr1-pr2-pr3-pr4 # n- shaped
  return (list(pr1=pr1,pr2=pr2,pr3=pr3,pr4=pr4,pr5=pr5,SHAPES = SHAPES))
}

