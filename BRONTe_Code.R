library(expm)
library(MASS)
library(truncnorm)
library(dplyr)
library(pracma)
library(abind)
# This function is for multi-way analysis (assuming rank R structure)

BRONTe <- function(X1,X2,XDiff,Y,N_sim=11000,burn_in=1000,Alpha,Beta,XTest1,XTest2,YTest,multi.source="yes",outcome="binary",zsd=1,rank=2){ #If multi.source=1 then multi-source analysis is performed, else, single-source
  
  # Progress bar 
  pb = txtProgressBar(min = 0, max = N_sim, initial = 0,style = 3) 
  N <- dim(X1)[1]        # Sample size
  Nn <- dim(XTest1)[1]   # Test set sample size
  P1 <- dim(X1)[2]        # Number of covariates
  P2 <- ifelse(length(dim(X1)) > 2,dim(X1)[3],1)        # Number of time points
  P3 <- dim(X1)[4]
  P <- dim(X2)[2]        # Number of covariates
  d <- ifelse(length(dim(X2)) > 2,dim(X2)[3],1)        # Number of time points
  R <- rank
  
  X2.1 <- X2[,,-XDiff]
  X2.2 <- X2[,,XDiff]

  y <- Y
  Xx1 <- XTest1
  Xx2 <- XTest2
  if(outcome=="binary"){
    N1  <- sum(y)         # Number of successes
    N0  <- N - N1         # Number of failures
  } else{
    N1  <- length(which(y>0))         # Number of successes
    N0  <- N - N1         # Number of failures
  }
  
  
  
  # Conjugate prior on the coefficients \beta ~ N(beta_0, Q_0)
  u1_0 <- rep(0, P1)      # u1 corresponds to covariates (p1)
  u2_0 <- rep(0, P2)      # u2 corresponds to time points (p2)
  u3_0 <- rep(0, P3)      # u3 corresponds to tissues (p3)
  
  # Initialize parameters
  tau2_v <- c(1)        # variance for v
  tau2_1 <- c(1)        # variance for w (source 1)
  tau2_a <- c(1)        # variance for w (source 2)
  y_var <- c(1)         # variance for y
  u1 <- rep(0, P1)        # initialize u1
  u2 <- rep(0, P2)        # initialize u2
  u3 <- rep(0, P3)        # initialize u3
  z  <- rep(0, length(y)) # initialize z
  
  # Matrices storing samples of the parameter
  u1_chain <- matrix(0, nrow = N_sim, ncol = P1*R)
  tauv_chain <- matrix(0, nrow = N_sim, ncol = 1)
  tau1_chain <- matrix(0, nrow = N_sim, ncol = 1)
  taua_chain <- matrix(0, nrow = N_sim, ncol = 1)
  u2_chain <- matrix(0, nrow = N_sim, ncol = P2*R)
  u3_chain <- matrix(0, nrow = N_sim, ncol = P3*R)
  prod_chain1 <- matrix(0, nrow = N_sim, ncol = (P1*P2*P3)) 
  y_var_chain <- matrix(0, nrow = N_sim, ncol = 1)
  
  #vecX <- matrix(X,nrow = N, ncol = P*d)
  U1s <- c(rep(0,P1*R))
  U1s_0 <- c(rep(0,P1*R))
  U2s <- c(rep(0,P2*R))
  U2s_0 <- c(rep(0,P2*R))
  U3s <- c(rep(0,P3*R))
  U3s_0 <- c(rep(0,P3*R))
  Xa <- matrix(X1,nrow = N, ncol = P1*R)
  Xb <- matrix(X1,nrow = N, ncol = P2*R)
  Xb1 <- matrix(X1,nrow = N, ncol = P2*R)
  Xc <- matrix(X1,nrow = N, ncol = P3*R)
  probtrain_chain <- matrix(0, nrow = N_sim, ncol = N) 
  probtest_chain <- matrix(0, nrow = N_sim, ncol = Nn) 
  probs_train <- Y
  probs_test <- c()
  U1s_r <- matrix(U1s,ncol = R)
  U2s_r <- matrix(U2s,ncol = R)  
  U3s_r <- matrix(U3s,ncol = R)
  
  # Conjugate prior on the coefficients \beta ~ N(beta_0, Q_0)
  w_0 <- rep(0, P*R)      # w corresponds to covariates (p)
  v_0 <- rep(0, d*R)      # v corresponds to time points (d)
  
  # Initialize parameters
  tau2_v <- c(1)        # variance for v
  tau2_1 <- c(1)        # variance for w (source 1)
  tau2_2 <- c(1)        # variance for w (source 2)
  y_var <- c(1)         # variance for y
  w <- rep(0, P*R)        # initialize w
  v <- rep(0, d*R)        # initialize v
  
  # Matrices storing samples of the parameter
  w_chain <- matrix(0, nrow = N_sim, ncol = P*R)
  tauv_chain <- matrix(0, nrow = N_sim, ncol = 1)
  tau1_chain <- matrix(0, nrow = N_sim, ncol = 1)
  tau2_chain <- matrix(0, nrow = N_sim, ncol = 1)
  v_chain <- matrix(0, nrow = N_sim, ncol = d*R)
  prod_chain2 <- matrix(0, nrow = N_sim, ncol = (P*d)) 
  y_var_chain <- matrix(0, nrow = N_sim, ncol = 1)
  
  vecX <- matrix(X2,nrow = N, ncol = P*d)
  Ws <- c(rep(0,P*R))
  Ws_0 <- c(rep(0,P*R))
  Vs <- c(rep(0,d*R))
  Vs_0 <- c(rep(0,1*R))
  U2s.2 <- c(rep(0,1*R))
  Xa2 <- matrix(X2.1,nrow = N, ncol = P*R)
  Xb2 <- matrix(X2.1,nrow = N, ncol = P2*R)
  probtrain_chain <- matrix(0, nrow = N_sim, ncol = N) 
  probtest_chain <- matrix(0, nrow = N_sim, ncol = Nn) 
  probs_train <- Y
  probs_test <- c()
  X1_est <- c()
  X2_est <- c()
  Ws_r <- matrix(Ws,ncol = R)
  Vs_r <- matrix(Vs,ncol = R)
  est_B1_list <- list()
  est_B2_list <- list()
  rbp <- matrix(nrow = P2,ncol = R)
  
  if(multi.source=="yes"){  
    for (t in 2:N_sim) {
      
      
      
      mat.fun.1 <- function(x,dim1=P1,dim2=(P2*P3)){matrix(c(x),nrow = dim1,ncol = dim2)}
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.1(X1[i,,,])
      }
      
      U2U3_list <- matrix(nrow = R,ncol = P2*P3)
      for(r in 1:R){
        U2U3_list[r,] <- as.vector(matrix(U2s,nrow = P2)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      
      
      for(i in 1:N){
        Xa[i,] <- as.vector(Xn[[i]] %*% t(U2U3_list))
      }
      
      
      # Update Mean of z
      mu_z <- Xa %*% U1s 
      # Draw latent variable z from its full conditional: z | W, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of W
      Q_0 <- diag(c(rep(tau2_1, P1*R)))
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P1*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa, Xa)))
      
      # Alt method to help with convergence
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa, Xa))))
      M1 <- Var %*% (prec_0 %*% U1s_0 + (yv)%*%crossprod(Xa, z)) 
      alt <- function(){
        c(mvrnorm(1, M1, Var_alt * s))
      }
      U1s <- tryCatch( {c(mvrnorm(1, M1, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # Posterior variance, prior = IG(Alpha,Beta)
      tau2_1 <- 1/rgamma(1,Alpha + ((P1*R)/2) ,Beta + (1/2)*(sum((U1s)^2)))

      ############################################################     
      mat.list <- list()
      mat.fun.2 <- function(x,dim1=P2,dim2=P1,dim3=P3){
        size <- (P1*P2)
        tot <- length(c(x))
        num <- tot/size
        for(q in 1:num){
          a <- 1+(size*(q-1))
          b <- size*q
          mat.list[[q]] <- matrix(c(x)[a:b],nrow = dim1,ncol = dim2,byrow = TRUE)
        }
        do.call(cbind,mat.list)
      }
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.2(X1[i,,,])
      }
      
      
      U1U3_list <- matrix(nrow = R,ncol = P1*P3)
      for(r in 1:R){
        U1U3_list[r,] <- as.vector(matrix(U1s,nrow = P1)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      
      
      for(i in 1:N){
        Xb1[i,] <- as.vector(Xn[[i]] %*% t(U1U3_list))
      }
      
      
      
      for(a in 1:N){
        for(r in 1:R){
          rbp[,r] <-   Ws_r[,r] %*% X2.1[a,,]
        }
        Xb2[a,] <- as.vector(rbp)
      }

      Xb <- Xb1 + Xb2
      
      mu_z <- Xb %*% U2s 
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, P2*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P2*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb, Xb))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb, Xb))))
      # Compute posterior mean of V
      M2 <- Var %*% (prec_0 %*% U2s_0 + (yv)%*%crossprod(Xb, z)) 
      alt <- function(){
        c(mvrnorm(1, M2, Var_alt * s))
      }
      U2s <- tryCatch( {c(mvrnorm(1, M2, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      ########################################################3      
      
      mat.fun.3 <- function(x,dim1=P3,dim2=(P2*P1)){matrix(c(x),nrow = dim1,ncol = dim2,byrow = TRUE)}
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.3(X1[i,,,])
      }
      
      
      U1U2_list <- matrix(nrow = R,ncol = P1*P2)
      for(r in 1:R){
        U1U2_list[r,] <- as.vector(matrix(U1s,nrow = P1)[,r] %o% matrix(U2s,nrow = P2)[,r])
      }
      
      
      for(i in 1:N){
        Xc[i,] <- as.vector(Xn[[i]] %*% t(U1U2_list))
      }
      
      mu_z <- Xc %*% U3s 
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, P3*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P3*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xc, Xc))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xc, Xc))))
      # Compute posterior mean of V
      M3 <- Var %*% (prec_0 %*% U3s_0 + (yv)%*%crossprod(Xc, z)) 
      alt <- function(){
        c(mvrnorm(1, M3, Var_alt * s))
      }
      U3s <- tryCatch( {c(mvrnorm(1, M3, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      
      
      
      
      mbp <- matrix(nrow = N,ncol = R)
      for(r in 1:R){
        mbp[,r] <-  X2.2 %*%  Ws_r[,r]
      }
      Xb3 <- mbp
      
      
      
      mu_z <- Xb3 %*% U2s.2
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, 1*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,1*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb3, Xb3))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb3, Xb3))))
      # Compute posterior mean of V
      M <- Var %*% (prec_0 %*% Vs_0 + (yv)%*%crossprod(Xb3, z)) 
      alt <- function(){
        c(mvrnorm(1, M, Var_alt * s))
      }
      U2s.2 <- tryCatch( {c(mvrnorm(1, M, Var))}, #<---- Var * 10^-10
                         error = function(e){
                           alt()
                         }
      )
      Vs <- c(U2s.2, U2s)
      # v Variance fixed at 1
      tau2_v <- 1
      
      
      Vs_r <- matrix(Vs,ncol = R,byrow = TRUE)
      map <- matrix(nrow = P,ncol = R)
      
      for(c in 1:N){
        for(r in 1:R){
          map[,r] <-  X2[c,,] %*%  Vs_r[,r]
        }
        Xa2[c,] <- as.vector(map)
      }
      
      
      mu_z <- Xa2 %*% Ws
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_a, P*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa2, Xa2))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa2, Xa2))))
      # Compute posterior mean of V
      M <- Var %*% (prec_0 %*% Ws_0 + (yv)%*%crossprod(Xa2, z)) 
      alt <- function(){
        c(mvrnorm(1, M, Var_alt * s))
      }
      Ws <- tryCatch( {c(mvrnorm(1, M, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # Posterior variance, prior = IG(Alpha,Beta)
      tau2_a <- 1/rgamma(1,Alpha + ((P*R)/2) ,Beta + (1/2)*(sum((Ws)^2)))
      
      Ws_r <- matrix(Ws,ncol = R,byrow = TRUE)
      
      
      # Store the draws
      u1_chain[t, ] <- U1s
      tau1_chain[t, ] <- tau2_1
      u2_chain[t, ] <- U2s
      u3_chain[t, ] <- U3s
      for(r in 1:rank){
        est_B1_list[[r]] <- (matrix(U1s,nrow = P1)[,r] %o% matrix(U2s,nrow = P2)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      est_B1 <- Reduce('+', est_B1_list)
      prod_chain1[t, ] <- as.vector(est_B1)
      
      
      b_param <- sum((Y - probs_train)^2)
      y_var <- ifelse(zsd=1,1,1/rgamma(1,0.001 + N/2, 0.001 + b_param))
      
      y_var_chain[t,] <- y_var
      # Store the draws
      w_chain[t, ] <- Ws
      taua_chain[t, ] <- tau2_a
      tau1_chain[t, ] <- tau2_1
      v_chain[t, ] <- Vs
      for(r in 1:rank){
        est_B2_list[[r]] <- (matrix(Ws,nrow = P)[,r] %o% matrix(Vs,nrow = d)[,r])
      }
      est_B2 <- Reduce('+', est_B2_list)
      prod_chain2[t, ] <- as.vector(est_B2)
      for(n in 1:N){
        probs_train[n] <- if(outcome=="binary"){pnorm((X1[n,,,] %*% est_B1) + sum(dot(X2[n,,], est_B2)) ) }else{ (X1[n,,,] %*% est_B1) + sum(dot(X2[n,,] %*% est_B2))  } 
      }
      for(n in 1:Nn){
        probs_test[n] <- if(outcome=="binary"){pnorm((Xx1[n,,,] %*% est_B1) + sum(dot(Xx2[n,,], est_B2)) ) }else{ (Xx1[n,,,] %*% est_B1) + sum(dot(Xx2[n,,] %*% est_B2))  }
      }
      probtrain_chain[t, ] <- probs_train
      probtest_chain[t, ] <- probs_test
      
      
      y_var_chain[t,] <- y_var
      setTxtProgressBar(pb,t)
    }
  }else{
    for (t in 2:N_sim){
      
      
      
      mat.fun.1 <- function(x,dim1=P1,dim2=(P2*P3)){matrix(c(x),nrow = dim1,ncol = dim2)}
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.1(X1[i,,,])
      }
      
      U2U3_list <- matrix(nrow = R,ncol = P2*P3)
      for(r in 1:R){
        U2U3_list[r,] <- as.vector(matrix(U2s,nrow = P2)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      
      
      for(i in 1:N){
        Xa[i,] <- as.vector(Xn[[i]] %*% t(U2U3_list))
      }
      
      
      # Update Mean of z
      mu_z <- Xa %*% U1s 
      # Draw latent variable z from its full conditional: z | W, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of W
      Q_0 <- diag(c(rep(tau2_1, P1*R)))
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P1*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa, Xa)))
      
      # Alt method to help with convergence
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa, Xa))))
      M1 <- Var %*% (prec_0 %*% U1s_0 + (yv)%*%crossprod(Xa, z)) 
      alt <- function(){
        c(mvrnorm(1, M1, Var_alt * s))
      }
      U1s <- tryCatch( {c(mvrnorm(1, M1, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # Posterior variance, prior = IG(Alpha,Beta)
      tau2_1 <- 1/rgamma(1,Alpha + (((P+P1)*R)/2) ,Beta + (1/2)*(sum(c(U1s,Ws)^2)))
      ############################################################     
      mat.list <- list()
      mat.fun.2 <- function(x,dim1=P2,dim2=P1,dim3=P3){
        size <- (P1*P2)
        tot <- length(c(x))
        num <- tot/size
        for(q in 1:num){
          a <- 1+(size*(q-1))
          b <- size*q
          mat.list[[q]] <- matrix(c(x)[a:b],nrow = dim1,ncol = dim2,byrow = TRUE)
        }
        do.call(cbind,mat.list)
      }
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.2(X1[i,,,])
      }
      
      
      U1U3_list <- matrix(nrow = R,ncol = P1*P3)
      for(r in 1:R){
        U1U3_list[r,] <- as.vector(matrix(U1s,nrow = P1)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      
      
      for(i in 1:N){
        Xb1[i,] <- as.vector(Xn[[i]] %*% t(U1U3_list))
      }
      
      
      
      for(a in 1:N){
        for(r in 1:R){
          rbp[,r] <-   Ws_r[,r] %*% X2.1[a,,]
        }
        Xb2[a,] <- as.vector(rbp)
      }

      Xb <- Xb1 + Xb2
      
      mu_z <- Xb %*% U2s 
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, P2*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P2*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb, Xb))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb, Xb))))
      # Compute posterior mean of V
      M2 <- Var %*% (prec_0 %*% U2s_0 + (yv)%*%crossprod(Xb, z)) 
      alt <- function(){
        c(mvrnorm(1, M2, Var_alt * s))
      }
      U2s <- tryCatch( {c(mvrnorm(1, M2, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      ########################################################3      
      
      mat.fun.3 <- function(x,dim1=P3,dim2=(P2*P1)){matrix(c(x),nrow = dim1,ncol = dim2,byrow = TRUE)}
      Xn <- list()
      for(i in 1:N){
        Xn[[i]] <- mat.fun.3(X1[i,,,])
      }
      
      
      U1U2_list <- matrix(nrow = R,ncol = P1*P2)
      for(r in 1:R){
        U1U2_list[r,] <- as.vector(matrix(U1s,nrow = P1)[,r] %o% matrix(U2s,nrow = P2)[,r])
      }
      
      
      for(i in 1:N){
        Xc[i,] <- as.vector(Xn[[i]] %*% t(U1U2_list))
      }
      
      mu_z <- Xc %*% U3s 
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, P3*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P3*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xc, Xc))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xc, Xc))))
      # Compute posterior mean of V
      M3 <- Var %*% (prec_0 %*% U3s_0 + (yv)%*%crossprod(Xc, z)) 
      alt <- function(){
        c(mvrnorm(1, M3, Var_alt * s))
      }
      U3s <- tryCatch( {c(mvrnorm(1, M3, Var))}, #<---- Var * 10^-10
                       error = function(e){
                         alt()
                       }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      
      
      
      
      Ws_r <- matrix(Ws,ncol = R)
      
      
      
      mbp <- matrix(nrow = N,ncol = 2)
      for(r in 1:R){
        mbp[,r] <-  X2.2 %*%  Ws_r[,r]
      }
      Xb3 <- mbp
      
      
      
      mu_z <- Xb3 %*% U2s.2
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, 1*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,1*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb3, Xb3))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb3, Xb3))))
      # Compute posterior mean of V
      M <- Var %*% (prec_0 %*% Vs_0 + (yv)%*%crossprod(Xb3, z)) 
      alt <- function(){
        c(mvrnorm(1, M, Var_alt * s))
      }
      U2s.2 <- tryCatch( {c(mvrnorm(1, M, Var))}, #<---- Var * 10^-10
                         error = function(e){
                           alt()
                         }
      )
      Vs <- c(U2s.2, U2s)
      # v Variance fixed at 1
      tau2_v <- 1
      
      
      Vs_r <- matrix(Vs,ncol = R)
      map <- matrix(nrow = P,ncol = R)
      
      for(c in 1:N){
        for(r in 1:R){
          map[,r] <-  X2[c,,] %*%  Vs_r[,r]
        }
        Xa2[c,] <- as.vector(map)
      }
      
      
      mu_z <- Xa2 %*% Ws
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, P*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,P*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa2, Xa2))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa2, Xa2))))
      # Compute posterior mean of V
      M <- Var %*% (prec_0 %*% Ws_0 + (yv)%*%crossprod(Xa2, z)) 
      alt <- function(){
        c(mvrnorm(1, M, Var_alt * s))
      }
      Ws <- tryCatch( {c(mvrnorm(1, M, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1

      
      
      
      # Store the draws
      u1_chain[t, ] <- U1s
      tau1_chain[t, ] <- tau2_1
      u2_chain[t, ] <- U2s
      u3_chain[t, ] <- U3s
      for(r in 1:rank){
        est_B1_list[[r]] <- (matrix(U1s,nrow = P1)[,r] %o% matrix(U2s,nrow = P2)[,r] %o% matrix(U3s,nrow = P3)[,r])
      }
      est_B1 <- Reduce('+', est_B1_list)
      prod_chain1[t, ] <- as.vector(est_B1)
      
      b_param <- sum((Y - probs_train)^2)
      y_var <- ifelse(zsd=1,1,1/rgamma(1,0.001 + N/2, 0.001 + b_param))
      
      y_var_chain[t,] <- y_var
      # Store the draws
      w_chain[t, ] <- Ws
      taua_chain[t, ] <- tau2_a
      tau1_chain[t, ] <- tau2_1
      v_chain[t, ] <- Vs
      for(r in 1:rank){
        est_B2_list[[r]] <- (matrix(Ws,nrow = P)[,r] %o% matrix(Vs,nrow = d)[,r])
      }
      est_B2 <- Reduce('+', est_B2_list)
      prod_chain2[t, ] <- as.vector(est_B2)
      for(n in 1:N){
        probs_train[n] <- if(outcome=="binary"){pnorm((X1[n,,,] %*% est_B1) + sum(dot(X2[n,,], est_B2)) ) }else{ (X1[n,,,] %*% est_B1) + sum(dot(X2[n,,] %*% est_B2))  } 
      }
      for(n in 1:Nn){
        probs_test[n] <- if(outcome=="binary"){pnorm((Xx1[n,,,] %*% est_B1) + sum(dot(Xx2[n,,], est_B2)) ) }else{ (Xx1[n,,,] %*% est_B1) + sum(dot(Xx2[n,,] %*% est_B2))  }
      }
      probtrain_chain[t, ] <- probs_train
      probtest_chain[t, ] <- probs_test

      y_var_chain[t,] <- y_var
      setTxtProgressBar(pb,t)
    }
  }
  # Get posterior mean of v,w, product of v and w, test class probabilities
  post_w_msmw <- w_chain[-(1:burn_in), ]
  post_v_msmw <- v_chain[-(1:burn_in), ]
  post_prod <- prod_chain2[-(1:burn_in), ]
  post_probtest <- probtest_chain[-(1:burn_in), ]
  post_probtrain <- probtrain_chain[-(1:burn_in), ]
  # Get posterior mean of v,w, product of v and w, test class probabilities
  post_u1_msmw <- u1_chain[-(1:burn_in), ]
  post_u2_msmw <- u2_chain[-(1:burn_in), ]
  post_u3_msmw <- u3_chain[-(1:burn_in), ]
  post_prod_msmw <- prod_chain1[-(1:burn_in), ]
  post_tau1 <-   tau1_chain[-(1:burn_in), ]
  post_taua <-   taua_chain[-(1:burn_in), ]  
  
  
  
  #Correlations of coefficients with truth and misclassification rate
  if(outcome=="binary"){
    class_tru <- YTest
    misclas <- class_tru - round(mean(post_probtest))
    misclassification_mse <- sum(abs(misclas))/Nn
  } else{
    class_tru <- YTest
    misclas <- class_tru - mean(post_probtest)
    misclassification_mse <- sum((misclas)^2)/sum((class_tru^2))
  }
  
  if(outcome=="binary"){
    class_tru <- Y
    misclas <- class_tru - round(colMeans(post_probtrain))
    train_mse <- sum(abs(misclas))/Nn
  } else{
    class_tru <- Y
    misclas <- class_tru - colMeans(post_probtrain)
    train_mse <- sum((misclas)^2)/sum((class_tru^2))
  }
  
  
  output <- list(post_w_msmw,post_v_msmw,post_u1_msmw,post_u2_msmw,post_u3_msmw,post_prod,post_prod_msmw,misclassification_mse,train_mse,post_probtest,post_tau1,post_taua)
  names(output) <- c("Ws","Vs","U1s","U2s","U3s","Source 2 prod","Source 1 prod","Misclass/MSE","TrainErr","PosteriorTest","Tau1","TauA")
  return(output)  
}
