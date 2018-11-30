library(gtools)
library(coda)
library(MASS)
library(invgamma)


density_estimation <- function(data,m,N.iter,n.iter,N.burn) {
  
  n=nrow(data)
  p=ncol(data)
  
  # initialize parameters
  K <- matrix(0,N.iter,n)
  theta <- matrix(0,N.iter,m)
  alpha <- matrix(0,N.iter,p)
  beta <- matrix(0,N.iter,p)
  
  theta[1,]<-rep(1/m,m)
  beta[1,]<-rep(1.5,p)
  alpha[1,]<-rep(1.5,p)
  K[1,] <- sample.int(m,size = n, replace = TRUE)
  
  
  # subset to get a dataset whose all dimensions are observed
  part_data=data[complete.cases(data), ]
  part_n=nrow(part_data)
  a = rep(part_n**0.4+1, p); b = apply(part_data, 2, var, na.rm=TRUE)
  
  # compute the knots 
  knots = matrix(0,p,m)
  x=part_data[,1]
  index_order=round((1:(m-2))*part_n/(m-1))
  index_order=ifelse(index_order==0,1,index_order)
  knots[1,] = c(min(x),sort(x)[index_order],max(x));
  for (i in 2:p) {
    y <- part_data[,i]
    y.ord=y[order(x)]
    index_order=round((1:(m-2))*part_n/(m-1))
    index_order=ifelse(index_order==0,1,index_order)
    knots[i,]=y.ord[c(1,index_order,part_n)]
  }
  
  
  
  # If no missing values, only use Gibbs_sampler to sample parameters
  if(part_n==n){
    
    list=Gibbs_sampler(data, m, alpha[1,], beta[1,], K[1,], theta[1,], a,b,knots, N.iter, N.burn)
    
    return(list(data=data,theta=list$theta,alpha=list$alpha,beta=list$beta,K=list$K,knots=knots))
  }
  
  
  # If the data has missing values, use missing_value_sampler and Gibbs_sampler
  if(part_n!=n){
    
    # record the location where data is missing
    miss_ind = matrix(0, n, p)
    for (i in 1:n){
      miss_ind[i,] <- is.na(data[i,])
    }
    
    wholedata=matrix(0, n, p)
    
    for(l in 1:(N.iter-1)){
      
      newdata=missing_value_sampler(data, miss_ind, knots, theta[l,], alpha[l,], beta[l,])
      wholedata=wholedata+newdata
      para=Gibbs_sampler(newdata, m, alpha[l,], beta[l,], K[l,], theta[l,], a,b,knots, n.iter, n.iter-1)
      
      alpha[l+1,]=para$alpha
      beta[l+1,]=para$beta
      theta[l+1,]=para$theta
      K[l+1,]=para$K
      
    }  
    
    return(list(newdata=wholedata/(N.iter-1),theta=theta[(N.burn+1):N.iter,],alpha=alpha[(N.burn+1):N.iter,],beta=beta[(N.burn+1):N.iter,],K=K[(N.burn+1):N.iter,],knots=knots))
    
  }
}


# density_estimation_old

density_estimation_old <- function(data,m,N.iter,n.iter,N.burn) {
  
  n=nrow(data)
  p=ncol(data)
  
  # initialize parameters
  K <- matrix(0,N.iter,n)
  theta <- matrix(0,N.iter,m)
  alpha <- matrix(0,N.iter,p)
  
  theta[1,]<-rep(1/m,m)
  alpha[1,]<-rep(1.5,p)
  K[1,] <- sample.int(m,size = n, replace = TRUE)
  
  
  # subset to get a dataset whose all dimensions are observed
  part_data=data[complete.cases(data), ]
  part_n=nrow(part_data)
  a = rep(part_n**0.4+1, p); b = apply(part_data, 2, var, na.rm=TRUE)
  #alpha[1,]<-(apply(data, 2, var, na.rm=TRUE)/n**(2/5))^0.5
  # compute the knots 
  knots = matrix(0,p,m)
  x=part_data[,1]
  index_order=round((1:(m-2))*part_n/(m-1))
  index_order=ifelse(index_order==0,1,index_order)
  knots[1,] = c(min(x),sort(x)[index_order],max(x));
  for (i in 2:p) {
    y <- part_data[,i]
    y.ord=y[order(x)]
    index_order=round((1:(m-2))*part_n/(m-1))
    index_order=ifelse(index_order==0,1,index_order)
    knots[i,]=y.ord[c(1,index_order,part_n)]
  }
  
  
  
  # If no missing values, only use Gibbs_sampler to sample parameters
  if(part_n==n){
    
    list=fix_beta_2(data, m, alpha[1,], K[1,], theta[1,], a,b, knots, N.iter, N.burn)
    
    return(list(data=data,theta=list$theta,alpha=list$alpha,beta=list$beta,K=list$K,knots=knots))
  }
  
  
  # If the data has missing values, use missing_value_sampler and Gibbs_sampler
  if(part_n!=n){
    
    # record the location where data is missing
    miss_ind = matrix(0, n, p)
    for (i in 1:n){
      miss_ind[i,] <- is.na(data[i,])
    }
    
    wholedata=matrix(0, n, p)
    
    for(l in 1:(N.iter-1)){
      
      newdata=missing_value_sampler_lastyear(data, miss_ind,knots, theta[l,], alpha[l,])
      wholedata=wholedata+newdata
      para=fix_beta_2(newdata, m, alpha[l,], K[l,], theta[l,],a,b,knots, n.iter, n.iter-1)
      
      
      alpha[l+1,]=para$alpha
      theta[l+1,]=para$theta
      K[l+1,]=para$K
      
    }  
    
    return(list(newdata=wholedata/(N.iter-1),theta=theta[(N.burn+1):N.iter,],alpha=alpha[(N.burn+1):N.iter,],K=K[(N.burn+1):N.iter,],knots=knots))
    
  }
}


