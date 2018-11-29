library(coda)
library(invgamma)
library(gtools)
library(foreach)
library(doSNOW)
library(doParallel)
set.seed(740)
# A hierarchial model to generate 2-dimensional data
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
data1 = cbind(x,y)

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

n.cl <- detectCores() # dectect the number of the cores in the local PC
cl <- makeCluster(n.cl-1) # use (n.cl-1) cores in the machine

cross_validation=function(data, nfolder,cores){
  # sd for each colunm
  # std = apply(data, 2, sd, na.rm=TRUE)
  ## nfolder: the number of folders
  #the length of true data
  n=dim(data)[1]
  #the dim of parameter
  p=dim(data)[2]

  size1=floor(n/nfolder)
  #size of the complete traning data
  size_train=n-size1


  index=1:n
  index1=index
  #split the complete data into nfolder group
  group_index=list(nfolder)
  for(i in 1:(nfolder-1)){
    group_index[[i]]=sample(index1,size1,replace=FALSE)
    index1=index1[!index1 %in% group_index[[i]]]
  }
  group_index[[nfolder]]=index1

  squared_error=function(index,m){
    test_data1=data[index,]
    if(size1==1)
      test_data1=matrix(test_data1,ncol=p)
    train_data=data[-index,]

    b=Gibbs_sampler(train_data,m,N.burn=5, N.iter=20)
    knots=b$knots
    alpha = b$alpha
    beta = b$beta
    theta = b$theta
    K = b$K
    #the number of test data
    n1=nrow(test_data1)
    sum1=0
    index2=1:p

    for(i in 1:n1){
      for(j in 1:p){
        real_value = test_data1[i,j]
        data1 = data[i,-j]
        if(p == 2){
          knots1 = matrix(knots[-j,],1)
        }
        else(
          knots1 = knots[-j,]
        )
        predict_value = c()
        for(l in 1: dim(beta)[1]){
          alpha1 = alpha[l,-j]
          beta1 = beta[l,-j]
          prod = 1
          weight = rep(0,m)
          weight1 = rep(0,m)
          for(k in 1:m){
            prod = 1
            for(s in 1: length(alpha1)){
              prod = prod * beta1[s]/(2 * alpha1[s]*gamma(1/beta1[s])) * exp(-abs((data1[s]-knots1[s,k])/alpha1[s])^beta1[s])
            }
            weight[k] = prod * theta[k]
            weight1[k] = weight[k]/sum(weight)
          }
          predict_value=c(predict_value, sum(weight1 * knots[j,]))
          #      print(l)
        }
        predict_value = mean(predict_value)
        sum1 = sum1 + (predict_value - real_value)^2/var(data[,j])
      }
    }
    return(sum1)


    # for(i in 1:n1){
    #   index_complete=index2[!is.na(test_data1[i,])]
    #   real_value=test_data1[i,index_complete]
    #   if(length(index_complete)==1){
    #     predict_value=apply(theta,1,function(x) sum(x*knots[index_complete,]))
    #     predict_value=mean(predict_value)
    #     sum1=sum1+((predict_value-real_value)/std[index_complete])^2
    #   }
    #   else{
    #     for(j in 1:length(index_complete)){
    #       k=index_complete[j]
    #       real_value=test_data1[i,k]
    #       if(length(index_complete)==2){
    #         knots1=matrix(knots[index_complete[-j],],1)
    #         predict_value = c()
    #         for(l in 1:dim(lambdas)[1]) {
    #           lambdas1=lambdas[l, index_complete[-j]]
    #           weight1=apply(knots1,2,function(x) prod(dnorm(test_data1[i,index_complete[-j]],mean=x,sd=lambdas1)))
    #           predict_value=c(predict_value, sum(weight1*theta[l,]/sum(weight1*theta[l,])*knots[k,]))
    #         }
    #       }
    #       else{
    #         knots1=knots[index_complete[-j],]
    #         predict_value = c()
    #         for(l in 1:dim(lambdas)[1]) {
    #           lambdas1=lambdas[l, index_complete[-j]]
    #           weight1=apply(knots1,2,function(x) prod(dnorm(test_data1[i,index_complete[-j]],mean=x,sd=lambdas1)))
    #           predict_value=c(predict_value, sum(weight1*theta[l,]/sum(weight1*theta[l,])*knots[k,]))
    #         }
    #       }
    #       predict_value=mean(predict_value)
    #       sum1=sum1+((predict_value-real_value)/std[index_complete[j]])^2
    #


    # if(length(index_complete)==2){
    #   knots1=matrix(knots[index_complete[-j],],1)
    #
    #   lambdas1=lambdas[index_complete[-j]]
    # }
    # else{
    #   knots1=knots[index_complete[-j],]
    #
    #   lambdas1=lambdas[index_complete[-j]]
    # }
    # weight1=apply(knots1,2,function(x) prod(dnorm(test_data1[i,index_complete[-j]],mean=x,sd=lambdas1)))
    # predict_value=apply(theta,1,function(x) sum(weight1*x/sum(weight1*x)*knots[k,]))
    # predict_value=mean(predict_value)
    # sum1=sum1+(predict_value-real_value)^2
  }

  #choose the possible value of m
  center=ceiling(size_train/log(size_train))
  seq1=(center-5):(center+5)
  step=floor((size_train-center-2)/3)
  seq3=seq(center+2+step,center+2+step*3,step)
  step=floor((center-3-2)/4)
  seq2=seq(center-2-4*step,center-2-step,step)
  m1=unique(c(3,4,seq2,seq1,seq3))
  m_num = length(m1)
  print(m1)
  Sys.sleep(10)

  registerDoSNOW(cores)
  errors <- foreach(i = 1:m_num,
                    .combine = 'cbind') %:%
    foreach(j = 1:nfolder, .combine = c) %dopar% {
      return(squared_error(group_index[[j]], m1[i]))
    }

  errors2 <- apply(errors,2,mean)
  return(m1[which.min(errors2)])
}
#################################### Parallel Computing ####################################

# #t1 <- Sys.time()
# registerDoSNOW(cl)
# errors <- foreach(i = 1:m_num,
#                   .combine = 'cbind',
#                   .export=c("density_estimate_new",
#                             "conditional_sampler",
#                             "posterior_Gibbs_sampler",
#                             "rdirichlet",
#                             "posterior_Gibbs_sampler.lambda_sq")) %:%
#   foreach(j = 1:nfolder, .combine = c) %do% {
#     return(squared_error(group_index[[j]], m1[i]))
#   }
# #t2 <- Sys.time()
# #print(difftime(t2,t1))

############################################################################################

#   err <- rep(0,m_num)
#   for (i in 1:m_num){
#     for (j in 1:nfolder){
#       err[i] = err[i] + squared_error(group_index[[j]],m1[i])
#     }
#     print(i)
#   }
#   return(m1[which.min(err)])
# }

# t1 <- Sys.time()
# cv_miss_new3(aa[[2]], 5, cl)
# t2 <- Sys.time()
# print("new")
# print(difftime(t2, t1)) # Time difference of 27.13687 mins
#
# t1 <- Sys.time()
# cv_miss_new(aa[[2]], 5)
# t2 <- Sys.time()
# print("old")
# print(difftime(t2, t1)) # Time difference of 2.331492 hours
cross_validation(data1,5,cl)
