#' elnet_coord Function
#' @description 
#' 
#' This function  fits elastic net to data using coordinate descent algorithm. The simulation setting is specified in homework 1 problem 3
#' part 2.
#' @usage 
#' 
#' elnet_coord()
#' @value The function will return a numeric vector of regression coefficients and together with a solution path plot.
#' @export


elnet_coord=function(){
  library("MASS")
  # Data generation
  p=20
  n=20
  alpha=0
  lambda_grid=seq(from=3,to=0,length.out = 100)
  beta_0=matrix(c(2,0,-2,0,1,0,-1,c(rep(0,13))), ncol = 1)
  sig=diag(c(rep(1,p)))
  sig[1,2]=0.8
  sig[2,1]=0.8
  sig[5,6]=0.8
  sig[6,5]=0.8
  X=mvrnorm(n=n,mu=c(rep(0,p)),Sigma=sig)
  # data standardization
  for(l in 1: p){
    X[,l]=sqrt(n)*X[,l]/sqrt(sum(X[,l]^2))
  }
  
  E=as.matrix(rnorm(n,0,1))
  Y=X%*%beta_0+E
  
  
  
  
  # coordinate descent
  soft=function(z,gam){
    if(gam >= abs(z)){
      result=0
    }else if (z>0){
      result=z-gam
    } else{
      result=z+gam
    }
    return(result)
  }
  up=function(X,Y,beta,j,lambda){
    res=1/n*t(as.matrix(X[,j]))%*%(Y-X[,-j]%*%beta[-j])
    beta[j]=soft(res,lambda*alpha)/(1+lambda-lambda*alpha)
    return(beta)
  }
  codes=function(X,Y,lambda,alpha){
    dif=1
    beta_co=matrix(0,ncol = 1,nrow = p)
    while(dif>1e-4){
      beta_co_old=beta_co
      for(j in 1 : p){
        beta_co=up(X,Y,beta_co,j,lambda)
      }
      dif=sqrt(sum(beta_co-beta_co_old)^2)
    }
    return(beta_co)
  }
  
  
  
  
  beta_co=matrix(nrow = p, ncol=100)
  beta_co_l1=matrix(nrow = 1, ncol=100)
  beta_co_l0=matrix(nrow = 1, ncol=100)
  for(i in 1 : 100){
    beta_co[,i]=codes(X=X,Y=Y,lambda = lambda_grid[i],alpha=alpha)
    beta_co_l1[,i]=sum(abs(beta_co[,i]))
    beta_co_l0[,i]=sum(beta_co[,i]!=0)
  }
  
  matplot(x=t(beta_co_l1),y=t(beta_co),type = "l",lty=1,xlab = "L1 Norm",ylab = "Coefficients",col = 1:20)
}