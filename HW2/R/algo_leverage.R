#' algo_leverage Function
#' @description 
#' 
#' This function implements algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows.
#' More specifically, the simulation setting of Figure 1 in https://onlinelibrary.wiley.com/doi/full/10.1002/wics.1324 is considered.
#' @usage 
#' 
#' algo_leverage()
#' @value a list contains two numeric vectors corresponding to uniform sampling (beta_unif) and leverage sampling (beta_lev). Each 
#' vector contains regression coefficients obtained using r = 10, 50, 100, 200, 300 in order.
#' @export

algo_leverage=function(){
  library("stats")
  X=as.matrix(rt(500,6,0))
  E=as.matrix(rnorm(500,0,1))
  Y=-X+E
  E_unif=matrix(nrow = 500,ncol = 5)
  E_lev=matrix(nrow = 500,ncol = 5)
  i=1
  for(r in c(10,50,100,200,300)){
    beta_unif=c()
    beta_lev=c()
    for(nos in 1:500){
      # Unif
      unif=sample(x=c(1:500),size = r,replace = F)
      X1=X[unif,]
      Y1=as.matrix(Y[unif,])
      beta_unif=c(beta_unif,solve(a=t(X1)%*%X1,b=t(X1)%*%Y1))
      
      # Leverage
      p_l=diag((X%*%t(X))/sum(X^2))
      s_5=sample(1:500,r,T,p_l)
      X_r5=X[s_5,]/sqrt(p_l[s_5])
      Y_r5=Y[s_5,]/sqrt(p_l[s_5])
      beta_lev=c(beta_lev,solve(a=t(X_r5)%*%X_r5,b=t(X_r5)%*%Y_r5))
    }
    i=i+1
  }
  return(list(beta_unif,beta_lev))
}