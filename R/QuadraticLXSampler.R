#.LXsample This function 
.LXsample_quadratic_2 <- function(t,y,JK,nsamps,nburn,anal_dtype,CONST,AXTREMA)#,NUM_CHG_PTS)
{
  #################################
  ##create the initial tree list / knot listp
  #######################0##########
  a1 = 0; a2 = 0; 
  if (AXTREMA == 1){
    a2 = -0.75; 
  }
  
  m.TREE = matrix(-1,nrow=6,3)
  m.TREE[1,1:3] = as.matrix(c(0,0.5,1))
  m.TREE[2,1:3] = as.matrix(c(0,0.25,0))
  m.TREE[3,1:3] = as.matrix(c(0,0.75,0))
  m.TREE[4,1:3] = as.matrix(c(0,0.25,0)) 
  m.TREE[5,1:3] = as.matrix(c(0, 1, 0))
  m.TREE[6,1:3] = as.matrix(c(0, 1, 0))
  CBX = shapespline2(m.TREE[1,]-.5, t-0.5,2)
  
  ONES = matrix(1,nrow=length(t),ncol=1)
  X = cbind(ONES,100*(CBX[,,3] - (a1+a2)*CBX[,,2]+ (a1*a2)*CBX[,,1]));
  BETAS = matrix(0,ncol(X),ncol=1)
  
  ptree = 0.5
  tau   = 0.25 
  x.Pred <- matrix(0,nrow=length(y),ncol=nsamps) #mean prediction of the data points
  lamT  = 1
  intLAM = 0.001
  p = 0.01
  
  omega <- matrix(1,nrow=nrow(y),ncol=1)
  r_disp = 1 ;
  LT <- rep(0,nrow(y))
  
  ###############################################################
  # Initialize the list for the simulated tempering
  ###############################################################
  BBETAS = vector("list",length(JK));
  BTAUS  = vector("list",length(JK)); 
  BA     = vector("list",length(JK)); 
  BX     = vector("list",length(JK));
  BCBX   = vector("list",length(JK)); 
  BTREE  = vector("list",length(JK)); 
  BCLAM  = vector("list",length(JK));
  BP     = vector("list",length(JK));
  BOMEGA = vector("list",length(JK));
  BZAUG  = vector("list",length(JK));
  BRDISP = vector("list",length(JK)); 
  ################################################################
  
  ################################################################
  #  print(omega)
  for (i in 1:length(JK)){
    BBETAS[[i]] = BETAS; 
    BTAUS[[i]]  = tau; 
    BCBX[[i]]   = CBX; 
    if (AXTREMA == 2){
      if (i == 2){
        a1 = -0.75;
        a2 =  0; 
      }
      if (i == 3){
        a1 = 0.75 ;
        a2 = 0; 
      }else{
        a1 = runif(1,-1,1)
        a2 = runif(1,-1,1)
      }
      
    }else{
      if (i == 2){
        a1 = 0.75 
      }
      if (i == 3){
        a1 = -0.75
      }
      if(i>3){
        a1 = runif(1,-1,1)
      }
    }
    
    BA[[i]] = c(a1, a2); 
    
    BX[[i]] = X; 
    BTREE[[i]] = m.TREE; 
    BCLAM[[i]] = lamT; 
    BP[[i]]    = p; 
    BOMEGA[[i]] = omega; 
    BRDISP[[i]] = r_disp; 
    BZAUG[[i]] =  switch(anal_dtype,y
                         ,0.5*(y-r_disp)*(1/omega)
                         ,(y[,1,drop=F]-0.5*y[,2,drop=F])*(1/omega))
  }
  
  MT  <- 1:(length(JK)-1)
  ADJ <- matrix(0,nrow=length(JK)-1,ncol=2)
  ADJ[,1] =  MT; ADJ[,2] =  MT+1; 
  ADJ = rbind(ADJ,c(1,3))
  ADJ = rbind(ADJ,c(1,4))
  ADJ = rbind(ADJ,c(2,4))
  ADJ = rbind(ADJ,c(2,5))
  
  
  rr = 1
  R_DISP   =  rep(0,nsamps)
  sa1      =  rep(0,nsamps)
  sa2      =  rep(0,nsamps)
  n_knots  =  rep(0,nsamps)
  h_lambda = rep(0,nsamps)
  h_p      = rep(0,nsamps)
  model_sample     = vector("list",nsamps);
  beta_sample      = vector("list",nsamps);
  knot_sample      = vector("list",nsamps); 
  tau_sample       = 1:nsamps
  
  
  ################################################################
  #
  #################################################################
  F.mat <- matrix(0,max(y)+1,max(y)+1) #matrix is only defined up to the maximum value of Y
  
  F.mat[1,1] <- 1
  
  for (m in 2:(max(y)+1)){
    for (j in 1:m){
      if ( j >= 2){
        F.mat[m,j] = F.mat[m-1,j]*(m-1)/m + (1/m)*F.mat[m-1,j-1]
      }else{
        F.mat[m,j] = F.mat[m-1,j]*(m-1)/m 
      }
    }
  }
  
  L <- rep(0,nrow(y))
  
  
  pb <- txtProgressBar(min=0, max=nsamps, style=3)
  
  
  for (i in 1:nsamps){
    
    setTxtProgressBar(pb,i)    
    
    if (runif(1) < 0.5){ # simple update 
      
      jj <- sample(1:length(JK),1)
      KK = JK[jj];
      BETAS = BBETAS[[jj]]; 
      tau = BTAUS[[jj]];
      CBX = BCBX[[jj]];  
      A = BA[[jj]]; a1 = A[1]; a2 = A[2]; 
      X = BX[[jj]]; 
      m.TREE = BTREE[[jj]]; 
      lamT = BCLAM[[jj]]; 
      
      
      #################################
      ##HYPER PRIOR LAMBDA
      betaT = BETAS[2:length(BETAS)]
      betaT = betaT[betaT > 0]
      
      ta = length(betaT) + .2
      tb = sum(betaT) + 2
      tg_ub = pgamma(1e-5,ta,tb)
      lamT = qgamma(runif(1,tg_ub,1),ta,tb)
      
      #other hyper prior on probability of a zero
      betaT = BETAS[2:length(BETAS)]
      V = betaT == 0
      t.n = length(betaT)
      ta = 2 + sum(V)
      tb = 18 + t.n - sum(V)
      p = rbeta(1,ta,tb)
      
      ####################################################################
      cbind(ONES,CONST*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      ####################################################################
      ###############################################################
      # Data Augmentation Step - 3 types of augmentations 
      #    (1) Normal -  nbinom
      #                 This is an identity augmentation and results in 
      #                 omega, the weight vector, being 1
      #    (2) Negative-Binomial -  'nbinom'
      #               - This is the augmentation based upon the paper 
      #               " Fully Bayesian Inference for neural models with-
      #                 negative binomial spiking" 2012 Nips
      #    (3) Binomial -   'bionom'
      #                   This is the augmentation based upon teh paper
      #                   "Bayesian infernce for logistic models using Polya-Gamma
      #                   latent variables." JASA 2013
      ###############################################################
      vecMean <- X%*%BETAS
      
      omega <- switch(anal_dtype,omega
                      ,.neg_binomial_augment(y,r_disp,vecMean)
                      ,.logig_binomial_augment(y,vecMean)); 
      
      augZ <- switch(anal_dtype,y
                     ,0.5*(y-r_disp)*(1/omega)
                     ,(y[,1,drop=F]-0.5*y[,2,drop=F])*(1/omega))
      
      
      if (anal_dtype ==2){
        ###############################################
        # Sample L
        ###############################################
        
        for( ii in 1:nrow(y)){
          if (y[ii] == 0){
            LT[ii] = 0;  
          }else{
            c <- seq(1,y[ii])
            temp <- r_disp^(1:y[ii])*F.mat[y[ii],c]
            
            LT[ii] <- as.numeric(sample(c,1, prob=temp/sum(temp)))
            
          }
        }
        
        ###############################################
        #
        #Sample r : To Do- Allow user to specify prior
        #
        ###############################################
        a = 5 + sum(LT)
        b = (5 - sum(X%*%BETAS - log((1+exp(X%*%BETAS)))))^(-1)
        r_disp <- rgamma(1,a,1/b)
        
        
        ###############################################  
      }else{
        r_disp = 1
      }
      
      if (anal_dtype == 1){ #only run if normal 
        tB <- 0.5*t(X%*%BETAS-augZ)%*%(X%*%BETAS-augZ)/KK+1
        tA <- 0.5*length(augZ)/KK+ 1
        tau <- as.vector((1/tB)*rgamma(1,tA,1))
      }else{
        tau <- 1
      }
      
      
      ###############################################################
      # GIBBS STEP 1: SAMPLE THE BETAS
      ################################################################
      BETAS = sampleBetas(augZ, X, BETAS,lamT, intLAM, p,tau/KK,omega)
      ################################################################ 
      ################################################################
      ################################################################
      # GIBBS STEP 2: SAMPLE THE CHANGE POINTS A1 AND A2
      # a1
      #    if(NUM_CHG_PTS == 2)
      ystar =  augZ - cbind(ONES,CONST*(CBX[,,3] - (a2)*CBX[,,2]))%*%BETAS;
      w = CONST*(-CBX[,,2]+(a2)*CBX[,,1])%*%BETAS[-1]
      a1  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5,omega)
      
      if (is.na(a1))
      { 
        print(jj)
        print(BETAS)
        #print(cbind(augZ,w,omega))
        
      }
      
      ################################################################
      # a2 
      ################################################################
      X = cbind(ONES,CONST*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      if (AXTREMA == 2){
        ystar =  augZ - cbind(ONES,CONST*(CBX[,,3] - (a1)*CBX[,,2]))%*%BETAS; 
        w = CONST*(-CBX[,,2]+(a1)*CBX[,,1])%*%BETAS[-1]
        a2  =  sampA(ystar,w,tau/KK,0,1,-0.5,0.5,omega)
      }  # a2 = 0.99
      
      # END STEP 2
      #################################################################
      X = cbind(ONES,CONST*(CBX[,,3] - (a1+a2)*CBX[,,2]+(a1*a2)*CBX[,,1]));
      #################################################################
      
      L = PROP.N.TREE(m.TREE);
      
      POS = L$pos;   
      
      if (L$INSERT==TRUE){
        ######################################################################################
        # GET THE NEW SHAPESPLINE
        NCBX = shapesplineInsertQuadratic(as.matrix(L$m.TREE[1,]-0.5,nrow=1), as.matrix(t,ncol=1)-0.5, 0.1,2,CBX,L$pos)
        NX = cbind(ONES,CONST*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        NBETAS = matrix(0,nrow=ncol(NX),ncol=1); 
        T_IDX <- (POS-1):(POS+2)+1;
        NBETAS[-T_IDX] = BETAS[-((POS-1):(POS+1)+1)]; 
        
        yStar = augZ - NX[,-T_IDX]%*%NBETAS[-T_IDX];
        LAM = lamT; 
        N_R = sinsertBeta_trivariate(as.matrix(yStar,ncol=1), NX[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(LAM,nrow=1),omega);
        
        ######################################################################################
        #Denominator Ratio for the equivalent delete
        ######################################################################################
        T_IDX <- (POS-1):(POS+1)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = 0*abs(m.TREE[1,]+10)*sqrt(((m.TREE[1,]-a1+10))^2*(m.TREE[1,]-a2+10)^2)+lamT;
       
          C_R =  sdeleteBeta_trivariate(yStar,X[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(LAM[T_IDX-1],ncol=1),omega)
        ##################################################################################
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB;
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        ##########################################################
        # Test statistic
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
        #########################################################
      }
      if (L$INSERT==FALSE){
     
        # GET THE NEW SHAPESPLINE
        NCBX = shapesplineDeleteQuadratic(as.matrix(L$m.TREE[1,]-0.5,ncol=1),as.matrix(t,ncol=1)-0.5,0.1,as.integer(2),CBX,as.integer(L$pos))
        NX = cbind(ONES,CONST*(NCBX[,,3] - (a1+a2)*NCBX[,,2]+(a1*a2)*NCBX[,,1]));
        
        #########################################################
        # Numerator of the Ration - Delete
        T_IDX <- (POS-1):(POS+1)+1;
        NBETAS = BETAS[-(POS+1),,drop=F];
        
        yStar = augZ - NX[,-T_IDX,drop=FALSE]%*%NBETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
         N_R = sdeleteBeta_trivariate(yStar,NX[,T_IDX,drop=FALSE],tau/KK,p,as.matrix(c(LAM,LAM,LAM),ncol=1),omega);
        #########################################################
        # Denomonator of the Ratio - Insert
        T_IDX <- (POS-1):(POS+2)+1;
        yStar = y - X[,-T_IDX,drop=FALSE]%*%BETAS[-T_IDX,,drop=FALSE];
        LAM = lamT
        
        C_R =  sinsertBeta_trivariate(yStar, X[,T_IDX,drop=FALSE], tau/KK, p,as.matrix(LAM),omega);
        #########################################################
        lnLike = N_R$LPROB; 
        ldLike = C_R$LPROB; 
        
        lnpMove = log(PROB.PROPOSAL(L$m.TREE,m.TREE))
        ldpMove = log(PROB.PROPOSAL(m.TREE,L$m.TREE))
        
        lnPrior = LOG.PROB.TREE(ptree,L$m.TREE)
        ldPrior = LOG.PROB.TREE(ptree, m.TREE)
        #########################################################
        # Test statistic
        test = lnLike + lnpMove + lnPrior - (ldLike + ldpMove + ldPrior); 
      }
      
      
      if ( test > 0 || runif(1) <exp(test))
      {
        m.TREE = L$m.TREE
        CBX   = NCBX; 
        X = NX;   
        
        # now sample the new beta
        if (L$INSERT){ 
          T_IDX <- (POS-1):(POS+2)+1;
          cmat = matrix(c( 0,0,0,0,   
                           1,0,0,0,   
                           0,2,0,0, 
                           0,0,3,0,
                           0,0,0,4,
                           1,2,0,0, 
                           1,0,3,0, 
                           1,0,0,4,
                           0,2,3,0,
                           0,2,0,4,
                           0,0,3,4,
                           0,2,3,4,
                           1,0,3,4,
                           1,2,0,4,
                           1,2,3,0,
                           1,2,3,4),ncol=4, nrow=16,byrow=T);
          tBETAS = c(0,0,0,0);
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS; 
          #  print(c(i,N_R$S_ELM)); 
        }else{
          T_IDX <- (POS-1):(POS+1)+1;  
          cmat = matrix(c( 0,0,0,   
                           1,0,0,   
                           0,2,0,
                           0,0,3,
                           1,2,0,
                           1,0,3,
                           0,2,3,
                           1,2,3) ,ncol=3, nrow=8,byrow=T);
          tBETAS = c(0,0,0); 
          tBETAS[cmat[N_R$S_ELM,]]= rtmvn(as.matrix(N_R$MEAN),as.matrix(N_R$VAR));
          NBETAS[T_IDX] = tBETAS;
          #  print(c(i,N_R$S_ELM)); 
        }          
        BETAS = NBETAS;     
      }
      #################################################################################
      BP[[jj]]  = p
      BBETAS[[jj]] = BETAS; 
      BTAUS[[jj]] = tau; 
      BCBX[[jj]] = CBX; 
      BA[[jj]] = c(a1,a2); 
      BX[[jj]] = X; 
      BTREE[[jj]] = m.TREE; 
      BCLAM[[jj]] = lamT; 
      BOMEGA[[jj]] = omega; 
      BZAUG[[jj]] = augZ; 
      BRDISP[[jj]] = r_disp; 
    }  else{
      #PERFORM A SWAP MOVE
      
      kk <- sample(1:nrow(ADJ),1)
      
      L1 = BCLAM[[ADJ[kk,1]]]
      L2 = BCLAM[[ADJ[kk,2]]]
      B1 = BBETAS[[ADJ[kk,1]]]
      B2 = BBETAS[[ADJ[kk,2]]]
      X1 = BX[[ADJ[kk,1]]]
      X2 = BX[[ADJ[kk,2]]]
      T1 = BTAUS[[ADJ[kk,1]]]
      T2 = BTAUS[[ADJ[kk,2]]]
      K1 = JK[[ADJ[kk,1]]]
      K2 = JK[[ADJ[kk,2]]]
      
      R1 = BRDISP[[ADJ[kk,1]]]
      R2 = BRDISP[[ADJ[kk,2]]]
      
      OMEG1 = BOMEGA[[ADJ[kk,1]]]
      OMEG2 = BOMEGA[[ADJ[kk,2]]]
      
      ZAUG1 <- BZAUG[[ADJ[kk,1]]]
      ZAUG2 <- BZAUG[[ADJ[kk,2]]]
      
      
      SS1 = (ZAUG1-X1%*%B1); SS1 = t(SS1*OMEG1)%*%SS1; 
      SS2 = (ZAUG2-X2%*%B2); SS2 = t(SS2*OMEG2)%*%SS2; 
      #note all of the prior probabilities are independent of K1 or K2 and thus
      #cancel in the ratio
      DEN  = -0.5*T1/K1*SS1  - 0.5*T2/K2*SS2 + 1/(K1*2)*sum(log(T1*(1/OMEG1)))+ 1/(K2*2)*sum(log(T2*(1/OMEG2))); 
      NUM  = -0.5*T1/K2*SS1  - 0.5*T2/K1*SS2 + 1/(K2*2)*sum(log(T1*(1/OMEG1)))+ 1/(K2*2)*sum(log(T2*(1/OMEG2)));#swap 
      ####################################################################
      test = NUM - DEN; 
      
      if (test >  0 || runif(1) < exp(test)){
        l = ADJ[kk,1]
        m = ADJ[kk,2]
        TB = BBETAS[[m]]; 
        TT = BTAUS[[m]];
        TC = BCBX[[m]]; 
        TA = BA[[m]];
        TX = BX[[m]];
        TTR= BTREE[[m]];
        TCL= BCLAM[[m]]; 
        TOMEGA = BOMEGA[[m]]
        TRDISP = BRDISP[[m]] 
        
        BBETAS[[m]]= BBETAS[[l]];
        BTAUS[[m]] = BTAUS[[l]];
        BCBX[[m]]  = BCBX[[l]]; 
        BA[[m]]    = BA[[l]];
        BX[[m]]    = BX[[l]];
        BTREE[[m]] = BTREE[[l]];
        BCLAM[[m]] = BCLAM[[l]];
        BOMEGA[[m]] = BOMEGA[[l]];
        BRDISP[[m]] = BRDISP[[l]];
        
        BBETAS[[l]] = TB; 
        BTAUS[[l]] = TT;
        BCBX[[l]] = TC; 
        BA[[l]] = TA;
        BX[[l]] = TX;
        BTREE[[l]] = TTR;
        BCLAM[[l]] = TCL; 
        BOMEGA[[l]] = TOMEGA; 
        BRDISP[[l]] = TRDISP; 
      }
    }  
    
    R_DISP[i]         = BRDISP[[rr]]
    tau_sample[i]     = BTAUS[[rr]]; 
    model_sample[[i]] = BTREE[[rr]]
    BETAS = BBETAS[[rr]]
    beta_sample[[i]] = BETAS
    h_p[i] = BP[[rr]]
    h_lambda[i] = BCLAM[[rr]]
    X = BX[[rr]]
    n_knots[i] = ncol(X)
    TEMP = X%*%BETAS; 
    
    t.mtree = BTREE[[rr]]; 
    beta_sample[[i]] = BETAS; 
    knot_sample[[i]] = t.mtree[1,]
    
    r = qcopy(x.Pred,TEMP,as.integer(nrow(y)),as.integer(i))
    TA = BA[[rr]]; 
    sa1[i] = TA[1];
    sa2[i] = TA[2];
  }
  
  
  # Estimate the mean 
  lmEST = rowMeans(x.Pred[,nburn:nsamps],na.rm=T)
  return(list(mcmcP = c(nburn,nsamps),lmEST=lmEST,sa1=sa1,sa2=sa2,
              beta_sample=beta_sample,model_sample=model_sample,
              x.Pred=x.Pred,dispersion=R_DISP,isIncreasing=(CONST > 0),EXTREMA = AXTREMA) ) 
}
