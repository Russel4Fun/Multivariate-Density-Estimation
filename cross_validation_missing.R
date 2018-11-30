library(foreach)
library(doSNOW)
library(doParallel)
n.cl <- detectCores() # dectect the number of the cores in the local PC
cl <- makeCluster(n.cl-1) # use (n.cl-1) cores in the machine
cross_validation_missing=function(data, nfolder, cl){
  # sd for each colunm
  std = apply(data, 2, sd, na.rm=TRUE)
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
  #split the complete data into 5 group
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

    b=density_estimation(train_data,m,N.iter = 3000, n.iter = 2000, N.burn = 1000)
    knots=b$knots
    theta=b$theta
    alpha = b$alpha
    beta = b$beta
    #the number of test data
    n1=nrow(test_data1)
    sum1=0
    index2=1:p
    for(i in 1:n1){
      index_complete=index2[!is.na(test_data1[i,])]
      real_value=test_data1[i,index_complete]
      if(length(index_complete)==1){
        predict_value=apply(theta,1,function(x) sum(x*knots[index_complete,]))
        predict_value=mean(predict_value)
        sum1=sum1+((predict_value-real_value)/std[index_complete])^2
      }
      else{
        for(j in 1:length(index_complete)){
          data1 = test_data1[i,index_complete[-j]]
          k=index_complete[j]
          real_value=test_data1[i,k]
          if(length(index_complete)==2){
            knots1=matrix(knots[index_complete[-j],],1)
            predict_value = c()
            for(l in 1: dim(beta)[1]){
              alpha1 = alpha[l,index_complete[-j]]
              beta1 = beta[l,index_complete[-j]]
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
          }
          else{
            knots1 = knots[index_complete[-j],]
            predict_value = c()
            for(l in 1: dim(beta)[1]){
              alpha1 = alpha[l,index_complete[-j]]
              beta1 = beta[l,index_complete[-j]]
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
            }
          }
          predict_value=mean(predict_value)
          sum1=sum1+((predict_value-real_value)/std[index_complete[j]])^2
        }
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

  #t1 <- Sys.time()
  registerDoSNOW(cl)
  errors <- foreach(i = 1:m_num,
                    .combine = 'cbind',
                    .export=c("density_estimation",
                              "missing_value_sampler",
                              "rinvgamma",
                              "rdirichlet",
                              "Gibbs_sampler")) %:%
    foreach(j = 1:nfolder, .combine = c) %dopar% {
      return(squared_error(group_index[[j]], m1[i]))
    }
  #t2 <- Sys.time()
  #print(difftime(t2,t1))
  for(i in 1:length(errors)){
    print(errors[i])
  }
  error2 <- apply(errors, 2, mean)
  return(m1[which.min(error2)])
}


cross_validation_missing_fixbeta2 <-function(data,nfolder,cl){
  # sd for each colunm
  std = apply(data, 2, sd, na.rm=TRUE)
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
  #split the complete data into 5 group
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

    b=density_estimation_old(train_data,m,N.iter = 3000, n.iter = 2000, N.burn = 1000)
    knots=b$knots
    theta=b$theta
    alpha = b$alpha
    #the number of test data
    n1=nrow(test_data1)
    sum1=0
    index2=1:p
    for(i in 1:n1){
      index_complete=index2[!is.na(test_data1[i,])]
      real_value=test_data1[i,index_complete]
      if(length(index_complete)==1){
        predict_value=apply(theta,1,function(x) sum(x*knots[index_complete,]))
        predict_value=mean(predict_value)
        sum1=sum1+((predict_value-real_value)/std[index_complete])^2
      }
      else{
        for(j in 1:length(index_complete)){
          data1 = test_data1[i,index_complete[-j]]
          k=index_complete[j]
          real_value=test_data1[i,k]
          if(length(index_complete)==2){
            knots1=matrix(knots[index_complete[-j],],1)
            predict_value = c()
            for(l in 1: dim(alpha)[1]){
              alpha1 = alpha[l,index_complete[-j]]
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
          }
          else{
            knots1 = knots[index_complete[-j],]
            predict_value = c()
            for(l in 1: dim(alpha)[1]){
              alpha1 = alpha[l,index_complete[-j]]
              prod = 1
              weight = rep(0,m)
              weight1 = rep(0,m)
              for(k in 1:m){
                prod = 1
                for(s in 1: length(alpha1)){
                  prod = prod/(alpha1[s]*gamma(1/2)) * exp(-abs((data1[s]-knots1[s,k])/alpha1[s])^2)
                }
                weight[k] = prod * theta[k]
                weight1[k] = weight[k]/sum(weight)
              }
            }
          }
          predict_value=mean(predict_value)
          sum1=sum1+((predict_value-real_value)/std[index_complete[j]])^2
        }
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

  #t1 <- Sys.time()
  registerDoSNOW(cl)
  errors <- foreach(i = 1:m_num,
                    .combine = 'cbind',
                    .export=c("density_estimation_old",
                              "missing_value_sampler_lastyear",
                              "rinvgamma",
                              "rdirichlet",
                              "fix_beta_2")) %:%
    foreach(j = 1:nfolder, .combine = c) %dopar% {
      return(squared_error(group_index[[j]], m1[i]))
    }
  #t2 <- Sys.time()
  #print(difftime(t2,t1))
  for(i in 1:length(errors)){
    print(errors[i])
  }
  error2 <- apply(errors, 2, mean)
  return(m1[which.min(error2)])
}
