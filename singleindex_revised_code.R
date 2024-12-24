library(quantreg)
library(fda)
library(rqPen)
library(tmvtnorm)
library(KernSmooth)
set.seed(1)
tau <- 0.5; c <- 10^{-6}

#generate data
p1 <- 100; p2 <- 100; n <- 200
mu <- rep(0,p1+p2)
sigma <- matrix(0,p1+p2,p1+p2)
for(i in 1:(p1+p2)){
  for(j in 1:(p1+p2)){
    sigma[i,j] <- 0.5^abs(i-j)
  }
}
funf <- function(x){
  sum(abs(x))
}

x1g <- seq(-1,1,0.02)


mainf <- function(i){
  print(i)
  x1 <- rtmvnorm(n,mu,sigma,lower=rep(-2,p1+p2),upper=rep(2,p1+p2))
  #x1 <- rmvnorm(n,mu,sigma)
  x <- x1[,1:p1]; z <- x1[,(p1+1):(p1+p2)]
  beta1 <- c(1,1,1,rep(0,p1-3))/sqrt(3)
  beta2 <- c(2,1,rep(0,p2-2))/sqrt(5)
  max_iter <-500
  a <- rbinom(n,1,0.5)
  y <- c()
  for(i in 1:n){
    y[i] <- sin(x[i,]%*%beta1)*a[i]+exp(z[i,]%*%beta2)+0.5*rt(1,3)
    #y[i] <- (x[i,]%*%beta1)*(1-x[i,]%*%beta1)*a[i]+exp(z[i,]%*%beta2)+0.5*rt(1,3)
  }
  y <- y-qnorm(tau,sd=0.5)
  gg <- x*a
  init_beta <- coef(rq.pen(x=cbind(gg,z),y=y,tau=tau,penalty="LASSO",
                           lambda=c(0.04,0.1)))[-1,1]
  beta1h <- matrix(0,p1,max_iter); beta2h <- matrix(0,p2,max_iter)
  beta1h[,1] <- sort(abs(init_beta[1:p1]),decreasing=T)
  x <- x[,order(abs(init_beta[1:p1]),decreasing=T)]
  beta1h[,1] <- sign(beta1h[,1]) * beta1h[,1]/sqrt(sum(beta1h[,1]^2))
  beta2h[,1] <- sort(abs(init_beta[(p1+1):(p1+p2)]),decreasing=T)
  z <- z[,order(abs(init_beta[(p1+1):(p1+p2)]),decreasing=T)]
  beta2h[,1] <- sign(beta2h[,1]) * beta2h[,1]/sqrt(sum(beta2h[,1]^2))
  knots <- c(seq(1,n,by=ceiling(n/25)),n)
  delta <- matrix(0,max_iter,44)
  for(j in 2:max_iter){
    print(j)
    b1 <- splineDesign(x=x%*%beta1h[,j-1],knots=sort(x%*%beta1h[,j-1])[knots],ord=4, derivs=0,
                       outer.ok = T)*a
    b2   <- splineDesign(x=z%*%beta2h[,j-1],knots=sort(z%*%beta2h[,j-1])[knots],ord=4, derivs=0,
                         outer.ok = T)
    delta[j,] <- coef(rq(y~0+b1+b2,tau=tau))
    j1 <- rbind(beta1h[-1,j-1]/sqrt(1-sum(beta1h[-1,j-1]^2)),
                diag(p1-1))
    j2 <- rbind(beta2h[-1,j-1]/sqrt(1-sum(beta2h[-1,j-1]^2)),
                diag(p2-1))
    bd1 <- splineDesign(x=x%*%beta1h[,j-1],knots=sort(x%*%beta1h[,j-1])[knots],ord=4, derivs=1,
                        outer.ok = T)%*%delta[j,1:22]
    bd2 <- splineDesign(x=z%*%beta2h[,j-1],knots=sort(z%*%beta2h[,j-1])[knots],ord=4, derivs=1,
                        outer.ok = T)%*%delta[j ,23:44]
    val1 <- matrix(0,n,p1-1)
    val2 <- matrix(0,n,p2-1)
    for(i in 1:n){
      val1[i,] <- (bd1[i,]*x[i,]%*%j1)*a[i]
      val2[i,] <- bd2[i,]*z[i,]%*%j2
    }
    bc1 <- splineDesign(x=x%*%beta1h[,j-1],knots=sort(x%*%beta1h[,j-1])[knots],ord=4, derivs=0,
                        outer.ok = T)%*%delta[j,1:22]
    bc2 <- splineDesign(x=z%*%beta2h[,j-1],knots=sort(z%*%beta2h[,j-1])[knots],ord=4, derivs=0,
                        outer.ok = T)%*%delta[j ,23:44]
    int <- c()
    for(i in 1:n){
      int[i] <- a[i]*bc1[i]+bc2[i]-a[i]*bd1[i]*x[i,]%*%j1%*%beta1h[-1,j-1]-
        bd2[i]*z[i,]%*%j2%*%beta2h[-1,j-1]
    }
    beta_1<- coef(rq.pen.cv(x=cbind(val1,val2),y=(y-int),tau=0.5,penalty="SCAD"))
    #beta_1 <- coef(rq.pen(x=cbind(val1,val2),y=(y-int),tau=tau,
                          #penalty="LASSO",lambda=c(0.02,0.09)))[,2]
    if(sum(beta_1[2:p1]^2)<1){
      temp1 <-c(sqrt(1-sum(beta_1[2:p1]^2)),beta_1[2:p1])
    }
    else{
      beta_1[2:p1] <- (beta_1[2:p1]/sqrt(sum(beta_1[2:p1]^2)))*sqrt(1-c^2)
      temp1 <-c(sqrt(1-sum(beta_1[2:p1]^2)),beta_1[2:p1])
    }
    if(sum(beta_1[(p1+1):(p1+p2-1)]^2)<sqrt(1-c^2)){
      temp2 <- c(sqrt(1-sum(beta_1[(p1+1):(p1+p2-1)]^2)),
                 beta_1[(p1+1):(p1+p2-1)])
    }
    else{
      beta_1[(p1+1):(p1+p2-1)] <- (beta_1[(p1+1):(p1+p2-1)]/sqrt(sum(beta_1[(p1+1):(p1+p2-1)]^2)))*sqrt(1-c^2)
      temp2 <- c(sqrt(1-sum(beta_1[(p1+1):(p1+p2-1)]^2)),
                 beta_1[(p1+1):(p1+p2-1)])
    }
    if(sum(c(temp1,temp2)-c(beta1h[,j-1],beta2h[,j-1]))^2<0.00001){
      j <-j-1
      break
    }
    x <- x[,order(abs(temp1),decreasing=T)]
    z <- z[,order(abs(temp2),decreasing=T)]
    beta1h[,j] <- temp1[order(abs(temp1),decreasing=T)]
    beta2h[,j] <- temp2[order(abs(temp2),decreasing=T)]
  }
  tb1 <- ifelse(beta1==0,1,0)
  eb1 <- ifelse(beta1h[,j]==0,1,0)
  fp <- ifelse(beta1==0& beta1h[,j]!=0,1,0)
  fn <- ifelse(beta1!=0& beta1h[,j]==0,1,0)
  um <- uj <- rep(0,length(x1g))
  y_star <- y-splineDesign(x=z%*%beta2h[,j],knots=sort(z%*%beta2h[,j])[knots],ord=4, derivs=0,
                           outer.ok = T)%*%delta[j ,23:44]
  x_star <- x%*%beta1h[,j]
  #h_tau <- dpill(x_star,y_star)*(tau*(1-tau)/(dnorm(qnorm(0.5)^2)))^0.2
  h_tau <- 0.34
  g1 <- splineDesign(x=x%*%beta1h[,j],knots=sort(x%*%beta1h[,j])[knots],ord=4, derivs=0,
                     outer.ok = T)%*%delta[j,1:22]
  for(i in 1:length(x1g)){
    for(k in 1:nrow(x)){
      um[i] <- um[i]+dnorm((x1g[i]-as.vector(x[k,]%*%beta1h[,j]))/0.5,0,1)/h_tau
      uj[i] <- uj[i]+dnorm((x1g[i]-as.vector(x[k,]%*%beta1h[,j]))/0.5,0,1)/h_tau*
        dnorm((splineDesign(x=x1g[i],knots=sort(x%*%beta1h[,j-1])[knots],ord=4, derivs=0,
                            outer.ok = T)%*%delta[j,1:22]-g1[k])/0.5,0,1)/h_tau
    }
    um[i] <- um[i]/nrow(x)
    uj[i] <- uj[i]/nrow(x)
  }
  se <- (um/uj)*um^(-0.5)*h_tau^(-0.5)*n^(-0.5)*sqrt(tau*(1-tau))
  f <- function(z, h = h_tau) (3 /4*h) * (1 - (z / h)^2) * (abs(z) < h)
  K2 <- integrate(f,lower=-1,upper=1)
  fd <- function(x){(f(x,h=h_tau)*(-2*x/h_tau^2))^2}
  Kp2 <- integrate(fd,lower=-1,upper=1)
  ck <- Kp2$value/K2$value
  alpha_h <- sqrt(-2*log(h_tau))
  Q <- alpha_h+(1/alpha_h)*log(sqrt(ck)/2*pi)-log(-log(sqrt(1-0.05)))
  trt_effect_fun <- function(x0){
    kern <- c()
    for(i in 1:nrow(x)){
      kern[i] <- f((x0-as.vector(x[i,]%*%beta1h[,j])))/nrow(x)
    }
    optim_fun <- function(b){
      qloss <-sum(abs(y-b*a-splineDesign(x=z%*%beta2h[,j],knots=sort(z%*%beta2h[,j])[knots],ord=4, derivs=0,
                                         outer.ok = T)%*%delta[j ,23:44])*kern)
      return(qloss)
    }
    return(optimize(optim_fun,interval=c(-5,5)))
}
trt_effect_list <- sapply(x1g,trt_effect_fun)
trt_effect <- c()
for(i in 1:length(x1g)){
trt_effect[i] <- trt_effect_list[,i]$minimum
}
mse <- mean((as.vector(trt_effect)-x1g*(1-x1g))^2)
mae <- mean(abs(as.vector(trt_effect)-x1g*(1-x1g)))
#mse <- mean((as.vector(trt_effect)-sin(x1g))^2)
#mae <- mean(abs(as.vector(trt_effect)-sin(x1g)))
lb <- trt_effect-Q*se; ub <- trt_effect+Q*se
lbp <- trt_effect-qnorm(0.975)*se
ubp <- trt_effect+qnorm(0.975)*se
co <- ifelse(lb<= x1g*(1-x1g) & x1g*(1-x1g) <=ub,1,0)
cop <- ifelse(lbp <= x1g*(1-x1g) & x1g*(1-x1g) <=ubp,1,0)
return(list(trt_effect,mse,mae,mean(fp),mean(fn),lb,ub,co,cop))
}

rep <- 200
res <- lapply(1:rep,mainf)
mse <- mae <- fp <- fn <- com <- c()
for(i in 1:rep){
  mse[i] <- res[[i]][[2]]; mae[i] <- res[[i]][[3]]
  fp[i] <- res[[i]][[4]]; fn[i] <- res[[i]][[5]]
}
co <- cop <- matrix(0,rep,length(x1g))
for(i in 1:rep){
  co[i,] <- res[[i]][[8]]; cop[i,] <- res[[i]][[9]]
}
