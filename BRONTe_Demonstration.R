######################################
# Documentation
######################################

### Bayesian Regression On Numerous Tensors (BRONTe)

## Description

# Fits a Bayesian Linear model that can accommodate multi-source and multi-way data

## Usage

BRONTe(X.list,XDiff,Y.list,N_sim=11000,burn_in=1000,Alpha=1,Beta=1,X.Test.list,Y.Test.list,multi.source="yes",rank=2,outcome="binary")

## Arguments

# X.list- an object of class "list": data frames containing the independent variables in the model of dimension n x p x d where all data frames in the list are linked along the n and d dimensions.
# XDiff- A vector of values indicating the features unique to source 2 
# Y.list- an object of class "list": vectors containing the dependent variable of length n
# N_sim- number of iterations the Gibbs sampler will run
# burn_in- number of iterations to be removed as burn-in
# Alpha- shape parameter for inverse gamma hyperprior used on covariate effects, default is 1
# Beta- rate parameter for inverse gamma hyperprior used on covariate effects, default is 1
# X.Test.list- an object of class "list": data frames containing the test set of independent variables
# Y.Test.list- an object of class "list": vectors containing the dependent variable test set
# multi.source- logical; if TRUE, runs multi-source model
# rank- assumed rank of the covariate matrix. Must be a positive integer less than or equal to min(P3,d) 
# outcome- whether the dependent variable is "binary" or "continuous"

## Details

## Value
# returns a list with the following components:

#"Source 1 prod", the array of estimated coefficients from source 1
#"Source 2 prod", the array of estimated coefficients from source 2
#"TrainErr", measure of how the model fits the training data. Returns misclassification rate for binary data, and relative mean squared error for continuous data
#"Misclass/MSE", measure of how the model fits the test data. Returns misclassification rate for binary data, and relative mean squared error for continuous data
#"Tau1", variance estimates for the sources
#"TauA", variance estimates for the sources
#"U1s", chain of estimated way coefficients from source 1, used for prediction
#"U2s", chain of estimated way coefficients from source 1, used for prediction
#"U3s", chain of estimated way coefficients from source 1, used for prediction
#"Ws", chain of estimated way coefficients from source 2, used for prediction
#"Vs", chain of estimated way coefficients from source 2, used for prediction


######################################
# Example
######################################




# Generate data X
N <- 100
P <- 10
P1 <- 5
P2 <- 5
d <- 2
errvar <- 1
x <- rnorm((N*(P)*(d)),0,1)
X <- array(c(x0), dim = c(N,P,d))

# True values of parameters
true_tau <- 1/rgamma(1,1,1)
true_tau1 <- 1/rgamma(1,1,1)
ratio <- .5

true_v <- rnorm(d,mean = 0, sd = sqrt(true_tau))  
true_w1 <- norm(P1,mean = 0, sd = sqrt(true_tau1))
true_w2 <- norm(P2,mean = 0, sd = ratio*sqrt(true_tau1))  

true_w <- c(true_w1,true_w2)

true_B <- as.vector(true_w %*% t(true_v))

mu <- c()
for(n in 1:N){
  mu[n] <- (t(true_w) %*% X[n,,] %*% t(t(true_v)))
}

# Generate observation data y
y <- rnorm(N,mean=mu,sd=rep(errvar,N))
Y <- y


# Generate Test data
xx0 <- rnorm((N*(P)*(d)),0,1)
Xtest <- array(c(xx0), dim = c(N,P,d))
mutest <- c()
for(n in 1:N){
  mutest[n] <- (t(true_w) %*% Xtest[n,,] %*% t(t(true_v)))
}
yy <- rnorm(N,mean=mu,sd=rep(errvar,N))
Ytest <- yy


source("BRONTe_Code_1.2.R") #See code file: https://github.com/BiostatsKim/BRONTe/blob/main/BRONTe_Code.R

Xlist1 <- X[,c(1:5),]
Xlist2 <- X[,c(6:10),]
Xlistfull <- list(Xlist1,Xlist2)

func.test <- BayesMSMW(X.list=Xlistfull,Y.list = Y,multi.source = "yes",rank = 1,outcome = "continuous")


func.test$TrainErr


func.pred <- predict.MSMW(Xtest=Xtest,R=1, Ws=func.test.2$Ws, Vs=func.test.2$Vs,outcome = "continuous")


misclas <- Ytest - func.pred$prediction
relative_mse <- sum((misclas)^2)/sum((class_tru^2))

