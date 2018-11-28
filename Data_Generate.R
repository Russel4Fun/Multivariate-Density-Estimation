set.seed(740)
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

