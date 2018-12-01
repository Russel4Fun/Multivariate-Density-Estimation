set.seed(740)
library(datasets)
library(mice)
library(MASS)
# A hierarchial model to generate 2-dimensional data
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
data1 = cbind(x,y)

# generating missing data
miss_x_prob = 0.2; miss_y_prob= 0.2
miss_x_index = sample(1:n, size = floor(miss_x_prob * n),replace = F)
miss_y_index = sample((1:n)[-miss_x_index], size = floor(miss_y_prob * n),replace = F)
x1 = x
x1[miss_x_index] = NA
y1 = y
y1[miss_y_index] = NA
missing_data1 <- cbind(x1,y1)


# another simulation data
miss_x_prob1 = 0.2; miss_y_prob1 = 0.2;miss_z_prob1 = 0.2
sigma <- matrix(data = c(1, 1.2, 2.2, 1.2, 1, 1.2, 2.2, 1.2, 1), nrow = 3)
sigma <- sigma%*%t(sigma)
data2 <- mvrnorm(n = 100, mu = c(5, 5, 5), Sigma = sigma)
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

# real windspeed data
data('iris')
data3 = as.matrix(iris[,1:4])
n3 = dim(data3)[1]
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
missing_data3 <- as.matrix(missing_data3)
