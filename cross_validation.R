# parallel processing
library(foreach)
library(doSNOW)
library(doParallel)
n.cl <- detectCores() # dectect the number of the cores in the local PC
cl <- makeCluster(n.cl-1) # use (n.cl-1) cores in the machine

cross_validation=function(data, nfolder){
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
        knots1 = knots[-j,]
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
          print(l)
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
  step=floor((center-3-2)/3)
  seq2=seq(center-2-3*step,center-2-step,step)
  m1=unique(c(seq2,seq1,seq3))
  m_num = length(m1)

  #t1 <- Sys.time()
  # registerDoSNOW(cl)
  # errors <- foreach(i = 1:m_num,
  #                   .combine = 'cbind') %:%
  #   foreach(j = 1:nfolder, .combine = c) %dopar% {
  #     return(squared_error(group_index[[j]], m1[i]))
  #   }
  #
  # return(errors)
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

  err <- rep(0,m_num)
  for (i in 1:m_num){
    for (j in 1:nfolder){
      err[i] = err[i] + squared_error(group_index[[j]],m1[i])
    }
    print(i)
  }
  return(m1[which.min(err)])
}

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
