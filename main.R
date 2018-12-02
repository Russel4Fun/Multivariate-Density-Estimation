#set.seed(740)
library(datasets)
library(mice)
library(MASS)
library(LaplacesDemon)

# A hierarchial model to generate 2-dimensional data
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
data1 = cbind(x,y)

# generating missing data 0.2
miss_x_prob = 0.2; miss_y_prob= 0.2
miss_x_index = sample(1:n, size = floor(miss_x_prob * n),replace = F)
miss_y_index = sample((1:n)[-miss_x_index], size = floor(miss_y_prob * n),replace = F)
x1 = x
x1[miss_x_index] = NA
y1 = y
y1[miss_y_index] = NA
data1_miss0.2 <- cbind(x1,y1)

# generating missing data 0.4
miss_x_prob = 0.4; miss_y_prob= 0.4
miss_x_index = sample(1:n, size = floor(miss_x_prob * n),replace = F)
miss_y_index = sample((1:n)[-miss_x_index], size = floor(miss_y_prob * n),replace = F)
x1 = x
x1[miss_x_index] = NA
y1 = y
y1[miss_y_index] = NA
data1_miss0.4 <- cbind(x1,y1)



###############################################################################################################
# another simulation data
sigma <- matrix(data = c(1,0.1,0.1,0.1,1,0.1,0.1,0.1,1),3,3)
sigma <- sigma%*%t(sigma)
mu1 = matrix(c(1,1,1),1,3)
data2 <- LaplacesDemon::rmvl(100,mu=mu1,Sigma = sigma)

# generating missing data 0.2
miss_x_prob1 = 0.2; miss_y_prob1 = 0.2;miss_z_prob1 = 0.2
missing_x_index_1 = sample.int(100, size = 100 * miss_x_prob1, replace = F)
missing_y_index_1 = sample.int(100, size = 100 * miss_y_prob1, replace = F)
missing_z_index_1 = sample.int(100, size = 100 * miss_z_prob1, replace = F)
missing_data2 <- data2
missing_data2[missing_x_index_1,1] = NA
missing_data2[missing_y_index_1,2] = NA
missing_data2[missing_z_index_1,3] = NA
bad = c()
for (i in 1:100){
  if(is.na(missing_data2[i,1]) & is.na(missing_data2[i,2]) & is.na(missing_data2[i,3])){
    bad = c(bad,i)
  }
}
if(length(bad) >0){
  missing_data2 = missing_data2[-bad,]
}
data2_miss0.2 = missing_data2

# generating missing data 0.4
miss_x_prob1 = 0.4; miss_y_prob1 = 0.4;miss_z_prob1 = 0.4
missing_x_index_1 = sample.int(100, size = 100 * miss_x_prob1, replace = F)
missing_y_index_1 = sample.int(100, size = 100 * miss_y_prob1, replace = F)
missing_z_index_1 = sample.int(100, size = 100 * miss_z_prob1, replace = F)
missing_data2 <- data2
missing_data2[missing_x_index_1,1] = NA
missing_data2[missing_y_index_1,2] = NA
missing_data2[missing_z_index_1,3] = NA
bad = c()
for (i in 1:100){
  if(is.na(missing_data2[i,1]) & is.na(missing_data2[i,2]) & is.na(missing_data2[i,3])){
    bad = c(bad,i)
  }
}
if(length(bad) >0){
  missing_data2 = missing_data2[-bad,]
}
data2_miss0.4 = missing_data2
###############################################################################################################
# real iris data
data('iris')
data3 = as.matrix(iris[,1:4])

n3 = dim(data3)[1]
# generating missing data 0.2
missing_a_prob = 0.2;missing_b_prob = 0.2;missing_c_prob = 0.2
missing_d_prob = 0.2
missing_a_index = sample.int(n3, size = n3 * missing_a_prob, replace = F)
missing_b_index = sample.int(n3, size = n3 * missing_b_prob, replace = F)
missing_c_index = sample.int(n3, size = n3 * missing_c_prob, replace = F)
missing_d_index = sample.int(n3, size = n3 * missing_d_prob, replace = F)
missing_data3 <- data3
missing_data3[missing_a_index,1] = NA
missing_data3[missing_b_index,2] = NA
missing_data3[missing_c_index,3] = NA
missing_data3[missing_d_index,4] = NA
bad1 = c()
for (i in 1:100){
  if(is.na(missing_data3[i,1]) & is.na(missing_data3[i,2]) & is.na(missing_data3[i,3]) & is.na(missing_data3[i,4])){
    bad1 = c(bad1,i)
  }
}
if(length(bad1) >0){
  missing_data3 = missing_data3[-bad1,]
}
data3_miss0.2 <- as.matrix(missing_data3)

# generating missing data 0.4
missing_a_prob = 0.4;missing_b_prob = 0.4;missing_c_prob = 0.4
missing_d_prob = 0.4
missing_a_index = sample.int(n3, size = n3 * missing_a_prob, replace = F)
missing_b_index = sample.int(n3, size = n3 * missing_b_prob, replace = F)
missing_c_index = sample.int(n3, size = n3 * missing_c_prob, replace = F)
missing_d_index = sample.int(n3, size = n3 * missing_d_prob, replace = F)
missing_data3 <- data3
missing_data3[missing_a_index,1] = NA
missing_data3[missing_b_index,2] = NA
missing_data3[missing_c_index,3] = NA
missing_data3[missing_d_index,4] = NA
bad1 = c()
for (i in 1:100){
  if(is.na(missing_data3[i,1]) & is.na(missing_data3[i,2]) & is.na(missing_data3[i,3]) & is.na(missing_data3[i,4])){
    bad1 = c(bad1,i)
  }
}
if(length(bad1) >0){
  missing_data3 = missing_data3[-bad1,]
}
data3_miss0.4 <- as.matrix(missing_data3)

###############################################################################################################
#data = data1
#data_miss0.2 = data1_miss0.2
#data_miss0.4 = data1_miss0.4
#data = data2
#data_miss0.2 = data2_miss0.2
#data_miss0.4 = data2_miss0.4
data = data3
data_miss0.2 = data3_miss0.2
data_miss0.4 = data3_miss0.4

n = dim(data)[1]
p = dim(data)[2]

m=6

###############################################################################################################
b <- density_estimation(data,m,N.iter=300,n.iter=300,N.burn = 200)
knots = b$knots
alpha = b$alpha
beta = b$beta
theta = b$theta

b_miss0.2 <- density_estimation(data_miss0.2,m,N.iter=300, n.iter=300,N.burn=200)
newdata_miss0.2 = b_miss0.2$newdata
knots_miss0.2 = b_miss0.2$knots
alpha_miss0.2 = b_miss0.2$alpha
beta_miss0.2 = b_miss0.2$beta
theta_miss0.2 = b_miss0.2$theta
for(i in 1:p){
  print(sum((data[,i]-newdata_miss0.2[,i])^2)/20)
}

b_miss0.4 <- density_estimation(data_miss0.4,m,N.iter=300,n.iter=300,N.burn=200)
newdata_miss0.4 = b_miss0.4$newdata
knots_miss0.4 = b_miss0.4$knots
alpha_miss0.4 = b_miss0.4$alpha
beta_miss0.4 = b_miss0.4$beta
theta_miss0.4 = b_miss0.4$theta
for(i in 1:p){
  print(sum((data[,i]-newdata_miss0.4[,i])^2)/40)
}
###############################################################################################################

b_old <- density_estimation_old(data,m,N.iter=300,n.iter=300,N.burn=200)
knots_old = b_old$knots
alpha_old = b_old$alpha
theta_old = b_old$theta


b_old_miss0.2 <- density_estimation_old(data_miss0.2,m,N.iter=300,n.iter=300,N.burn=200)
newdata_old_miss0.2 = b_old_miss0.2$newdata
knots_old_miss0.2 = b_old_miss0.2$knots
alpha_old_miss0.2 = b_old_miss0.2$alpha
theta_old_miss0.2 = b_old_miss0.2$theta
for(i in 1:p){
  print(sum((data[,i]-newdata_old_miss0.2[,i])^2)/20)
}

b_old_miss0.4 <- density_estimation_old(data_miss0.4,m,N.iter=300,n.iter=300,N.burn=200)
newdata_old_miss0.4 = b_old_miss0.4$newdata
knots_old_miss0.4 = b_old_miss0.4$knots
alpha_old_miss0.4 = b_old_miss0.4$alpha
theta_old_miss0.4 = b_old_miss0.4$theta
for(i in 1:p){
  print(sum((data[,i]-newdata_old_miss0.4[,i])^2)/40)
}
###############################################################################################################
fk <- function(x,knot,alpha,beta){
  y = beta/(2 * alpha*gamma(1/beta)) * exp(-abs((x-knot)/alpha)^beta)
  return(y)
}

###ks-test for random beta
alpha.estimate = apply(alpha, 2, mean)
beta.estimate = apply(beta, 2, mean)
theta.estimate = apply(theta, 2, mean)

ks_b = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate[k]*fk(x,knots[j,k],alpha.estimate[j],beta.estimate[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b[j])
}

###ks-test for random beta, missing r = 0.2
alpha.estimate_miss0.2 = apply(alpha_miss0.2, 2, mean)
beta.estimate_miss0.2 = apply(beta_miss0.2, 2, mean)
theta.estimate_miss0.2 = apply(theta_miss0.2, 2, mean)

ks_b_miss0.2 = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate_miss0.2[k]*fk(x,knots_miss0.2[j,k],alpha.estimate_miss0.2[j],beta.estimate_miss0.2[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b_miss0.2[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b_miss0.2[j])
}

###ks-test for random beta, missing r = 0.4
alpha.estimate_miss0.4 = apply(alpha_miss0.4, 2, mean)
beta.estimate_miss0.4 = apply(beta_miss0.4, 2, mean)
theta.estimate_miss0.4 = apply(theta_miss0.4, 2, mean)

ks_b_miss0.4 = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate_miss0.4[k]*fk(x,knots_miss0.4[j,k],alpha.estimate_miss0.4[j],beta.estimate_miss0.4[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b_miss0.4[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b_miss0.4[j])
}

###############################################################################################################
fk <- function(x,knot,alpha){
  y = 1/(alpha*gamma(1/2)) * exp(-abs((x-knot)/alpha)^2)
  return(y)
}
###ks-test for fix beta = 2
alpha.estimate_old = apply(alpha_old, 2, mean)
theta.estimate_old = apply(theta_old, 2, mean)

ks_b_old = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate_old[k]*fk(x,knots_old[j,k],alpha.estimate_old[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b_old[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b_old[j])
}

###ks-test for fix beta = 2, missing r = 0.2   
alpha.estimate_old_miss0.2 = apply(alpha_old_miss0.2, 2, mean)
theta.estimate_old_miss0.2 = apply(theta_old_miss0.2, 2, mean)

ks_b_old_miss0.2 = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate_old_miss0.2[k]*fk(x,knots_old_miss0.2[j,k],alpha.estimate_old_miss0.2[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b_old_miss0.2[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b_old_miss0.2[j])
}

###ks-test for fix beta = 2, missing r = 0.4
alpha.estimate_old_miss0.4 = apply(alpha_old_miss0.4, 2, mean)
theta.estimate_old_miss0.4 = apply(theta_old_miss0.4, 2, mean)

ks_b_old_miss0.4 = rep(0,p)
for(j in 1:p){
  xpdf<-function(x){
    sum=0
    for(k in 1:m){
      sum = sum + theta.estimate_old_miss0.4[k]*fk(x,knots_old_miss0.4[j,k],alpha.estimate_old_miss0.4[j])
    }
    return(sum)
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  ks_b_old_miss0.4[j]=ks.test(data[,j],"xcdf")$p
  print(ks_b_old_miss0.4[j])
}










######################################################################################################
##plot
fk <- function(x,knot,alpha,beta){
  y = beta/(2 * alpha*gamma(1/beta)) * exp(-abs((x-knot)/alpha)^beta)
  return(y)
}

single_den=function(x,theta,knots,alpha,beta){
  m = length(theta)
  sum = 0
  for(i in 1:m){
    sum = sum + theta[i]*fk(x,knots[i],alpha,beta)
  }
  return(sum)
}

summaries_new <- function(data, j, theta, knots, alpha, beta) {
  x=data[,j]
  x.grid=seq(min(x),max(x),l=length(x));
  fx.hat=matrix(0,length(x.grid),nrow(theta))
  for(i in 1:length(x.grid)){
    for(l in 1:nrow(theta)){
      fx.hat[i,l]=single_den(x.grid[i],theta[l,],knots[j,],alpha[l,j],beta[l,j])
    }
  }
  fx.hat.gibbs.l=apply(fx.hat,1,quantile,prob=0.025)
  fx.hat.gibbs=apply(fx.hat,1,mean)
  fx.hat.gibbs.u=apply(fx.hat,1,quantile,prob=0.975)
  return(list(fx.hat=fx.hat,fx.hat.gibbs.l=fx.hat.gibbs.l,fx.hat.gibbs=fx.hat.gibbs,fx.hat.gibbs.u=fx.hat.gibbs.u,x.grid=x.grid))
}

summaries_old <- function(data, j, theta, knots, alpha) {
  x=data[,j]
  x.grid=seq(min(x),max(x),l=length(x));
  fx.hat=matrix(0,length(x.grid),nrow(theta))
  for(i in 1:length(x.grid)){
    for(l in 1:nrow(theta)){
      fx.hat[i,l]=single_den(x.grid[i],theta[l,],knots[j,],alpha[l,j],2)
    }
  }
  fx.hat.gibbs.l=apply(fx.hat,1,quantile,prob=0.025)
  fx.hat.gibbs=apply(fx.hat,1,mean)
  fx.hat.gibbs.u=apply(fx.hat,1,quantile,prob=0.975)
  return(list(fx.hat=fx.hat,fx.hat.gibbs.l=fx.hat.gibbs.l,fx.hat.gibbs=fx.hat.gibbs,fx.hat.gibbs.u=fx.hat.gibbs.u,x.grid=x.grid))
}

par(mfrow=c(3,3))

for(i in 1:ncol(data)){
  s1=summaries_new(data, i, theta, knots, alpha, beta)
  s2=summaries_old(data, i, theta_old,knots_old,alpha_old)
  x.hmax=1.2*max(hist(data[,i],breaks=15,plot=F)$density)
  hist(data[,i],freq=F,breaks=15,ylim=c(0,x.hmax), xlab=NA,ylab= "Dnesity",main= NA,sub=paste0("(a",i,")"))
  lines(s1$x.grid,s1$fx.hat.gibbs,col="red",lwd=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.l,col="red",lty=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.u,col="red",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs,col="blue",lwd=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.l,col="blue",lty=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.u,col="blue",lty=2)
}

for(i in 1:ncol(data)){
  s1=summaries_new(newdata_miss0.2, i, theta_miss0.2, knots_miss0.2, alpha_miss0.2, beta_miss0.2)
  s2=summaries_old(newdata_old_miss0.2, i, theta_old_miss0.2,knots_old_miss0.2,alpha_old_miss0.2)
  x.hmax=1.2*max(hist(data[,i],breaks=15,plot=F)$density)
  hist(data[,i],freq=F,breaks=15,ylim=c(0,x.hmax), xlab=NA,ylab="Dnesity",main= NA,sub=paste0("(b",i,")"))
  lines(s1$x.grid,s1$fx.hat.gibbs,col="red",lwd=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.l,col="red",lty=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.u,col="red",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs,col="blue",lwd=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.l,col="blue",lty=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.u,col="blue",lty=2)
}

for(i in 1:ncol(data)){
  s1=summaries_new(newdata_miss0.4, i, theta_miss0.4, knots_miss0.4, alpha_miss0.4, beta_miss0.4)
  s2=summaries_old(newdata_old_miss0.4, i, theta_old_miss0.4,knots_old_miss0.4,alpha_old_miss0.4)
  x.hmax=1.2*max(hist(data[,i],breaks=15,plot=F)$density)
  hist(data[,i],freq=F,breaks=15,ylim=c(0,x.hmax), xlab=NA,ylab="Dnesity",main= NA,sub=paste0("(c",i,")"))
  lines(s1$x.grid,s1$fx.hat.gibbs,col="red",lwd=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.l,col="red",lty=2)
  #lines(s1$x.grid,s1$fx.hat.gibbs.u,col="red",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs,col="blue",lwd=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.l,col="blue",lty=2)
  #lines(s2$x.grid,s2$fx.hat.gibbs.u,col="blue",lty=2)
}

