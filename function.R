# Here the N.iter denotes the total number of iterations
library(invgamma)
library(gtools)
library(coda)


Gibbs_sampler <- function(data, m, N.iter, N.burn){
  n = dim(data)[1]
  p = dim(data)[2]
  a0 = rep(1,m)
  # initializing the MCMC
  K <- matrix(0,N.iter,n)
  theta <- matrix(0,N.iter,m)
  alpha <- matrix(0,N.iter,p)
  beta <- matrix(0,N.iter,p)
  theta[1,]<-rep(1/m,m)
  beta[1,]<-rep(1.5,p)
  alpha[1,]<-rep(1.5,p)
  K[1,] <- sample.int(m,size = n, replace = TRUE)
  knots = matrix(0,p,m)
  # matrix to restore the observed data values
  f = matrix(0,n,m)
  
  # choosing knots by sorting the data and setting equal length of step
  x = data[,1]
  index_order=round((1:(m-2))*n/(m-1))
  index_order=ifelse(index_order==0,1,index_order)
  knots[1,] = c(min(x),sort(x)[index_order],max(x));
  for (i in 2:p) {
    y <- data[,i]
    y.ord=y[order(x)]
    index_order=round((1:(m-2))*n/(m-1))
    index_order=ifelse(index_order==0,1,index_order)
    knots[i,]=y.ord[c(1,index_order,n)]
  }
  
  # Now we are doing the MCMC sampling
  for(l in 2:N.iter){
    for(k in 1:m){
      f[,k] = rep(1,n)
      for (j in 1:p){
        coeff <- 0.5 * beta[l-1,j]/(alpha[l-1,j] * gamma(1/beta[l-1,j]))
        pow <- -(abs((data[,j]-knots[j,k])/alpha[l-1,j])**beta[l-1,j])
        f[,k] <- f[,k] * coeff * exp(pow)
      }
    }
    
    theta.old <- theta[l-1,]
    # Now sample K
    for (i in 1:n){
      fvalue <- as.vector(f[i,])
      K[l,i] <- sample.int(m, size=1, prob = theta.old * fvalue, replace = TRUE)
    }
    
    # Now sample theta
    theta.new = rdirichlet(1,a0+tabulate(K[l,],nbins=m))
    theta[l,] <- theta.new
    
    # Now sample alpha with the help of inverse gamma distribution
    #a <- rep(n^0.4 + 1,p)
    # b <- rep(0,p)
    # for(j in 1:p){
    #   b[j] <- var(data[,j])
    # }
    a <- rep(2.01,p)
    b <- rep(1.01,p)
    
    alpha_to_beta <- rep(0,p)
    for (j in 1:p){
      shape = n/beta[l-1,j] + a[j]
      sum = 0
      for(i in 1:n){
        sum = sum + (abs(data[i,j] - knots[j,K[l,i]])) ** beta[l-1,j]
      }
      scale = sum + b[j]
      alpha_to_beta[j] <- rinvgamma(n = 1, shape = shape, rate = scale)
      alpha[l,j] = alpha_to_beta[j] ** (1/beta[l-1,j])
    }
    
    # finally we gonna sample beta, using rejection sampling
    for (j in 1:p){
      K0 <- function(beta_0){
        part1 <- (beta_0^(n+1))/(gamma(1/beta_0)^n) * (1/alpha[l-1,j])^(a[j]*beta_0)
        sum = 0
        for (i in 1:n){
          sum = sum + abs(data[i,j]-knots[j,K[l,i]]) ** beta_0
        }
        part2 <- exp(-(sum + b[j])/(alpha[l,j]^beta_0))
        result <- ifelse(1<=beta_0 & beta_0<=2, part1 * part2, 0)
        return(result)
      }
      
      K1 <- function(beta){
        coeff = (n+1)*log(beta)-n*log(gamma(1/beta))
        sum = 0
        for(i in 1:n){
          sum = sum + abs(data[i,j]-knots[j,K[l,i]]) ** beta
        }
        part2 = -(sum + b[j])/(alpha[l,j] ^ beta)
        result <- ifelse(1<= beta & beta<=2, exp(part2 + coeff), 0)
        return(result)
      }
      
      K2 <- function(beta){
        coeff <- (n+1) * log(beta) - n*log(gamma(1/beta)) - a[j] * beta * log(alpha[l,j])
        sum = 0
        for(i in 1:n){
          sum = sum + abs(data[i,j]-knots[j,K[l,i]]) ** beta
        }
        part2 = -(sum + b[j])/(alpha[l,j] ^ beta)
        result <- ifelse(1<= beta & beta<=2, part2 + coeff,0)
        return(result)
      }
      
      g <- function(beta){
        k <- alpha[l,j] ^ (-a[j])
        result <- ifelse(1<= beta & beta <= 2, k^beta, 0)
        return (result)
      }
      
      temp <- integrate(g,1,2)$value
      
      g1 <- function(beta){
        return(1/temp * g(beta))
      }
      
      g.inv <- function(beta){
        k <- alpha[l,j] ^ (-a[j])
        result <- log(temp * beta * log(k) +k)/log(k)
        return(ifelse(0<beta & beta<1, result, 0))
      }
      
      # M <- optimize(K1, interval = c(1,2), maximum = T)$objective
      M <- K1(1)
      
      Rej.samp <- function(max.try = 10000){
        try = 1
        while(try < max.try){
          u0 <- runif(1)
          # sample p from the g function
          samp <- g.inv(u0)
          u <- runif(1)
          if(log(u) + log(M) + log(g1(samp)) < K2(samp)){
            return(samp)
          }
          try = try + 1
        }
        return('increase max try')
      }
      beta[l,j] <- Rej.samp(100000)
    }
  }
  return(list(K = K[(N.burn+1):N.iter, ], theta = theta[(N.burn+1):N.iter, ], alpha = alpha[(N.burn+1): N.iter, ],
              beta = beta[(N.burn+1):N.iter, ], knots = knots))
}
