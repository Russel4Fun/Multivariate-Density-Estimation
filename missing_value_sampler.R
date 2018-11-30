set.seed(740)

#set true mean and sd functions
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}

#generate data:
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
full_data=cbind(x,y)

#random missing
miss_x_prop = 0.2; miss_y_prop= 0.2
miss_x_ind = sample(1:n, size =round(miss_x_prop*n),  replace = FALSE)
miss_y_ind = sample((1:n)[-miss_x_ind], size =round(miss_y_prop*n),  replace = FALSE)
x_miss = x; y_miss = y
x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA

#generate data that have missing value
miss_data=cbind(x_miss, y_miss)

#record the location where data is missing
miss_ind = matrix(0, n, p)
for (i in 1:n){
  miss_ind[i,] <- is.na(miss_data[i,])
}




#sample missing value given parameters

missing_value_sampler <- function(data, miss_ind, knots, theta, alpha, beta) {
  
  n=nrow(data)
  p=ncol(data)
  m=length(theta)
  
  complete_data=matrix(0,n,p)
  
  for (i in 1:n){
    
    miss_ind_i=c(1:p)[miss_ind[i,]==TRUE]
    nonmiss_ind_i=c(1:p)[miss_ind[i,]==FALSE]
    
    if(length(nonmiss_ind_i)==p){complete_data[i,]=data[i,]}
    
    else{
      complete_data[i,]=data[i,]
      
      freq=NULL
      for(k in 1:m){
        
        nomis=1
        for (j in nonmiss_ind_i){
          nomis=nomis*exp(-abs((data[i,j]-knots[j,k])/alpha[j])^beta[j])
        }
        freq=c(freq,theta[k]*nomis)
        
      }
      
      prob=freq/sum(freq)
      k=sample(1:m,1,prob=prob)
      
      for(j in miss_ind_i){
        　pdf=function(x){
            0.5*beta[j]*exp(-abs((x-knots[j,k])/alpha[j])^beta[j])/(alpha[j]*gamma(1/beta[j]))}
          cdf=function(x){
            sapply(x,function(x){integrate(pdf,lower=-Inf, upper=x)$value})}
          qf=function(u){
            h=function(x){cdf(x)-u} 
            uniroot(h,interval=c(knots[j,k]-10,knots[j,k]+10))$root}
          u=runif(1)
          while(u>0.99 | u<0.01){u=runif(1)}
          t=try(sapply(u,qf))
          if(isTRUE(class(t)=="try-error")) { complete_data[i,j]=mean(data, na.rm=TRUE) } 
          else {complete_data[i,j]=t} 
      }
    }  
  }
  
  return(complete_data)
  
}


# Last year method
missing_value_sampler_lastyear <- function(data, miss_ind, knots, theta, alpha) {
  
  n=nrow(data)
  p=ncol(data)
  m=length(theta)
  
  complete_data=matrix(0,n,p)
  
  for (i in 1:n){
    
    miss_ind_i=c(1:p)[miss_ind[i,]==TRUE]
    nonmiss_ind_i=c(1:p)[miss_ind[i,]==FALSE]
    
    if(length(nonmiss_ind_i)==p){complete_data[i,]=data[i,]}
    
    else{
      complete_data[i,]=data[i,]
      
      freq=NULL
      for(k in 1:m){
        
        nomis=1
        for (j in nonmiss_ind_i){
          nomis=nomis*dnorm(data[i,j],knots[j,k],alpha[j])
        }
        freq=c(freq,theta[k]*nomis)
        
      }
      
      prob=freq/sum(freq)
      k=sample(1:m,1,prob=prob)
      
      for(j in miss_ind_i){
        
       complete_data[i,j]=rnorm(1,knots[j,k],alpha[j]) 
      }
    }  
  }
  
  return(complete_data)
  
}
  
  