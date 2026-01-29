library(ggplot2)
library(tidyr)

##Parts a-d

N = 25
X = rnorm(mean = 0, sd = 1, n = N)

mean = c(0,0,2,2)
var = c(1/2, 2, 1/2, 2)


val = seq(-1,1, by = 0.01)
data = matrix(nrow = length(val), ncol = 5)

max = matrix(nrow = 1, ncol = 4)

colnames(data) = c("val", 1,2,3,4)
colnames(max) = c(1,2,3,4)
data[,1] = val

post = function(par){
  mu = par[1]
  
  l = -1/sqrt(2 * pi * tao_n)*exp(-(mu-mu_n)^2/(2*tao_n))
  
  return(l) 
  
}

for (i in 1:4){
  
  mu_n = (N/1*mean(X)+1/var[i]*mean[i])/(N/1+1/var[i])
  tao_n = 1/(N/1+1/var[i])
  
  f = 1/sqrt(2 * pi * tao_n)*exp(-(val-mu_n)^2/(2*tao_n))
  
  data[,i+1] = f
  
  ##assign(paste0("post_",i), f)
  
  maximum = optim(c(0), post)
  
  max[1,i] =  maximum$par
  
}

df = data.frame(data)
max

plot(df$val, df$X1,type = "l", col = "red")
lines(df$val, df$X2, type = "l", col = "blue")
lines(df$val, df$X3, type = "l", col = "green")
lines(df$val, df$X4, type = "l", col = "black")

legend(-1,2,legend=c("(0,1/2)","(0,2)","(2,1/2)", "(2,2)"), col=c("red","blue","green","black"),
       pch=c("l","l","l","l"))
