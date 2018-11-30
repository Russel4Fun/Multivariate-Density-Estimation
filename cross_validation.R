# parallel processing
library(foreach)
library(doSNOW)
library(doParallel)
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
    b = density_estimation(train_data,m,N.iter=300,n.iter=200,N.burn=100)
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
  }

  #choose the possible value of m
  seq1= 3:9
  seq2=seq(10,25,by=3)
  m1 = c(seq1,seq2)
  m_num = length(m1)
  print(m1)

  registerDoSNOW(cores)
  errors <- foreach(i = 1:m_num,
                    .combine = 'cbind',
                    .export = c('rdirichlet','rinvgamma','density_estimation','Gibbs_sampler','missing_value_sampler')) %:%
    foreach(j = 1:nfolder, .combine = c) %dopar% {
      return(squared_error(group_index[[j]], m1[i]))
    }

  errors2 <- apply(errors,2,mean)
  result <- m1[which.min(errors2)]
  return(result)
}

cross_validation_fixbeta2 = function(data, nfolder,cores){
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
    b= density_estimation_old(data = train_data, m = m, N.iter=300, n.iter = 200, N.burn=100)
    knots=b$knots
    alpha = b$alpha
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
        for(l in 1: dim(alpha)[1]){
          alpha1 = alpha[l,-j]
          prod = 1
          weight = rep(0,m)
          weight1 = rep(0,m)
          for(k in 1:m){
            prod = 1
            for(s in 1: length(alpha1)){
              prod = prod /(alpha1[s]*gamma(1/2)) * exp(-abs((data1[s]-knots1[s,k])/alpha1[s])^2)
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
  }

  seq1= 3:9
  seq2=seq(10,25,by=3)
  m1 = c(seq1,seq2)
  m_num = length(m1)
  print(m1)

  registerDoSNOW(cores)
  errors <- foreach(i = 1:m_num,
                    .combine = 'cbind',
                    .export = c('rdirichlet','rinvgamma','fix_beta_2',
                                'missing_value_sampler','density_estimation_old')) %:%
    foreach(j = 1:nfolder, .combine = c) %dopar% {
      return(squared_error(group_index[[j]], m1[i]))
    }

  errors2 <- apply(errors,2,mean)
  result <- m1[which.min(errors2)]
  return(result)
}

