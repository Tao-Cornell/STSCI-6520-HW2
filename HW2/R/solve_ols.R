#' solve_ols Function
#'
#' This function solves a linear system using Gauss-Seidel or Jacobi method.
#' @param X the coefficient matrix of the linear system and should be a square matrix.
#' @param b the column vector collects values of each linear equation. 
#' @param method the method is used to solve the linear system, either "Gauss" or "Jacobi" coresponding to Gauss-Seidel or Jacobi method respectively.
#' @param tol the convergence criterion which is measured by the \eqn{l_2} norm of the difference of between the current solution and the solution in the last iteration.
#' @return The solution vector of the linear system specified in the function  
#' @export


solve_ols=function(X,b,method,tol=1e-4){
  library(Matrix)
  L=tril(X)
  diag(L)=0
  U=triu(X)
  diag(U)=0
  D=Diagonal(x=diag(X))
  p=nrow(X)
  dif=tol+1
  iter_max=10^4
  iter=1
  x=matrix(0,nrow = p,ncol = 1)
  if(method == "Gauss"){
    while(dif>=tol && iter<=iter_max){
      x=solve(L+D)%*%(b-U%*%x)
      iter=iter+1
    } 
  }
  else if(method =="Jacobi"){
    while(dif>=tol && iter<=iter_max){
      x=solve(D)%*%(b-(L+U)%*%x)
      iter=iter+1
    } 
  }
  else{
    print("Method is not specified")
    x=NA
  }
  return(x)
}