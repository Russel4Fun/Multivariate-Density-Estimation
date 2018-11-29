set.seed(740)
library(datasets)
# A hierarchial model to generate 2-dimensional data
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
data1 = cbind(x,y)

# get the estimate using Gibbs Sampler
b <- Gibbs_sampler(data1,m=10,N.iter = 2000, N.burn = 1000)
dimension = dim(b$beta)[2]
beta.estimate = numeric(dimension)
alpha.estimate = numeric(dimension)
for (i in 1:dimension){
  beta.estimate[i] = mean((b$beta)[,i])
  alpha.estimate[i] = mean((b$alpha)[,i])
}

m.estimate <- cross_validation(data1,5)

# generating missing data
miss_x_prob = 0.4; miss_y_prob= 0.4
miss_x_index = sample.int(n, size = floor(miss_x_prob * n),replace = F)
miss_y_index = sample.int(n, size = floor(miss_y_prob * n),replace = F)
x1 = x
x1[miss_x_index] = NA
y1 = y
y1[miss_y_index] = NA
missing_data1 <- cbind(x1,y1)

# another simulation data
require(MASS)
miss_x_prob1 = 0.3; miss_y_prob1 = 0.3;miss_z_prob1 = 0.3
sigma <- matrix(data = c(1, 1.2, 2.2, 1.2, 1, 1.2, 2.2, 1.2, 1), nrow = 3)
sigma <- sigma%*%t(sigma)
data2 <- mvrnorm(n = 100, mu = c(5, 5, 5), Sigma = sigma)
missing_x_index_1 = sample.int(100, size = 100 * miss_x_prob1, replace = F)
missing_y_index_1 = sample.int(100, size = 100 * miss_y_prob1, replace = F)
missing_z_index_1 = sample.int(100, size = 100 * miss_z_prob1, replace = F)
missing_data2 <- data2
missing_data2[missing_x_index_1,1] = NA
missing_data2[missing_y_index_1,2] = NA
missing_data2[missing_x_index_1,3] = NA

# real iris data
data3 = iris[,1:4]
