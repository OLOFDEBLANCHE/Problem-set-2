library(ggplot2)
library(tidyr)

##Parts a-d

N = 100
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


##e-h

##Prior parameters 
mu = 0
tao = 1/2 
alpha = 1
beta = 1


##Initial values
mu_0= c(0,-2)
sigma_0 = c(1,5)

means = matrix(nrow = 2, ncol = 1001)
sigmas = matrix(nrow = 2, ncol = 1001)

for (j in 1:2){
  means[j,1] = mu_0[j]
  sigmas[j,1] = sigma_0[j]
  
  for(i in 1:1000){
    
    mu_n = (N/sigmas[j,i]*mean(X)+1/tao*mu)/(N/sigmas[j,i]+1/tao) ##Posterior mean
    tao_n = 1/(N/sigmas[j,i]+1/tao) ##Posterior variance
    
    means[j,i+1] = rnorm(n = 1, mean = mu_n, sd = sqrt(tao_n)) ##Draw next mean from distribution
    #conditional on data and variance[i]
    
    squares = (X - means[j,i+1])^2 
    
    shape_a = alpha + N/2 ##Posterior shape
    
    rate_b = beta + 1/2*sum(squares) ##Posterior rate
    
    sigmas[j,i+1] = 1/rgamma(n = 1,shape = shape_a, rate =rate_b) ##Draw next variance from
    #distribution conditional on data and mean[i + 1]
    
    
  } 
  
  

}

plot(1:1001, means[1,])
plot(1:1001, means[2,])
plot(1:1001, sigmas[1,])
plot(1:1001, sigmas[2,])

#f
dens_mean = density(means[1,])
dens_sigmas = density(sigmas[1,])


mu_n = (N/1*mean(X)+1/tao*mu)/(N/1+1/tao)
tao_n = 1/(N/1+1/tao)


f_mean = 1/sqrt(2 * pi * tao_n)*exp(-(dens_mean$x-mu_n)^2/(2*tao_n))

plot(dens_mean$x, dens_mean$y, type = "l")
lines(dens_mean$x, f_mean, type = "l", col = "blue")

plot(dens_sigmas$x, dens_sigmas$y, type = "l")



#h
dens2_mean = density(means[2,])
dens2_sigmas = density(sigmas[2,])


f_var = 1/sqrt(2 * pi * tao_n)*exp(-(dens2_mean$x[-1]-mu_n)^2/(2*tao_n))

plot(dens2_mean$x[-1], dens2_mean$y[-1], type = "l")
lines(dens2_mean$x[-1], f_var, type = "l", col = "blue")

plot(dens2_sigmas$x[-1], dens2_sigmas$y[-1], type = "l")




#i
##Repeated sampling of a VECTOR of the parameters. 
##Prior, joint distribution are assumed to be independent. 


##Initial values
mu_0 = 0
sigma_0 = 1

##Candidtate standard deviation
candidate_sd_mean = 1
candidate_sd_var = 1

means = rep(0, times = 1001)
vars = rep(0, times = 1001)

means[1] = mu_0
vars[1] = sigma_0


for(i in 1:1000){
  
  cand_mean = rnorm(n = 1, mean = means[i], sd = candidate_sd_mean) ##Draw mean from candidate dist
  cand_var = rnorm(n = 1, mean = vars[i], sd = candidate_sd_var) ##Draw var from candidate dist
  
  if(cand_var<0){
    cand_var = cand_var*(-1)
  }
  
  
  ##Candidate posterior
  mu_n_cand = (N/cand_var*mean(X)+1/tao*mu)/(N/cand_var+1/tao) ##Posterior mean of mu. 
  tao_n_cand = 1/(N/cand_var+1/tao) ##Posterior variance of mu
  
  f_cand_var = dgamma(cand_var, shape = alpha, rate = beta)
  
  f_cand_mean = dnorm(x = cand_mean, mean = mu_n_cand, sd = sqrt(tao_n_cand)) ##prob of cand mean
  #* likelihood of sample
  posterior_cand = f_cand_var * f_cand_mean
  

  ##Current posterior 
  mu_n = (N/vars[i]*mean(X)+1/tao*mu)/(N/vars[i]+1/tao) 
  tao_n = 1/(N/vars[i]+1/tao) 
  
  f_var = dgamma(vars[i], shape = alpha, rate = beta)

  f_mean = dnorm(x = means[i], mean = mu_n, sd = sqrt(tao_n)) ##Prob of current mean
  #* likelihood of sample 
  posterior_current = f_var * f_mean
  

  frac = posterior_cand / posterior_current
  
  
  u = runif(n = 1)
  
  if(u <= frac){
    means[i+1] = cand_mean
    vars[i+1] = cand_var
  }
  
  else{
    means[i + 1] = means[i]
    vars[i+1] = vars[i]
  }
  
}

##i
plot(1:1001, means)
plot(1:1001, vars)

dens3_mean = density(means)
dens3_sigmas = density(vars)


plot(dens3_mean$x, dens3_mean$y, type = "l")
lines(dens_mean$x, dens_mean$y, type = "l", col = "blue")

plot(dens_sigmas$x, dens_sigmas$y, type = "l")
lines(dens3_sigmas$x, dens3_sigmas$y, type = "l", col = "blue")


##j 

f_post = function(data, mu, sigma){
  
  probs = dnorm(data, mean=mu, sd = sqrt(sigma)) 
  f_mu = 1/4 * exp(-2*abs(mu))
  f_sigma = dexp(x = sigma, rate = 1)

  l = prod(probs)
  
  return(l*f_mu*f_sigma)
} 

##Initial values
mu_0 = 0
sigma_0 = 1

##Candidtate standard deviation
candidate_sd_mean = 1
candidate_sd_var = 1

means = rep(0, times = 1001)
vars = rep(0, times = 1001)

means[1] = mu_0
vars[1] = sigma_0


for(i in 1:1000){
  
  cand_mean = rnorm(n = 1, mean = means[i], sd = candidate_sd_mean) ##Draw mean from candidate dist
  cand_var = rnorm(n = 1, mean = vars[i], sd = candidate_sd_var) ##Draw var from candidate dist
  
  if(cand_var<0){
    cand_var = cand_var*(-1)
  }
  
  
  ##Candidate posterior

  posterior_cand = f_post(X, cand_mean, cand_var)
  
  
  ##Current posterior 

  posterior_current = f_post(X, means[i], vars[i])
  
  
  frac = posterior_cand / posterior_current
  
  
  u = runif(n = 1)
  
  if(u <= frac){
    means[i+1] = cand_mean
    vars[i+1] = cand_var
  }
  
  else{
    means[i + 1] = means[i]
    vars[i+1] = vars[i]
  }
  
}

##i
plot(1:1001, means)
plot(1:1001, vars)

dens4_mean = density(means)
dens4_sigmas = density(vars)

##Posterior density given known variance and data
mu_n = (N/1*mean(X)+1/tao*mu)/(N/1+1/tao)
tao_n = 1/(N/1+1/tao)
f_mean = 1/sqrt(2 * pi * tao_n)*exp(-(dens4_mean$x-mu_n)^2/(2*tao_n))
##

plot(dens4_mean$x, dens4_mean$y, type = "l")
lines(dens3_mean$x, dens3_mean$y, type = "l", col = "red")
lines(dens_mean$x, dens_mean$y, type = "l", col = "blue")
lines(dens4_mean$x, f_mean, type = "l", col = "green")
legend(-0.1,8,legend=c("MH-pr2","MH-pr1","Gibbs", "Known var"), col=c("black","red", "blue","green"),
       pch=c("l","l","l","l"))


plot(dens_sigmas$x, dens_sigmas$y, type = "l")
lines(dens4_sigmas$x, dens3_sigmas$y, type = "l", col = "blue")










