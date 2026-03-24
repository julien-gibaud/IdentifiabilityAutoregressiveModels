library(MARSS)

simulationAR1 <- function(n=1000, T=50, rho, s2eps, s2eta, x0){
  list.rho <- c()
  list.s2eps <- c()
  list.s2eta <- c()
  pb <- txtProgressBar(char=' \u27a4', style=3)
  for(i in 1:n){
    set.seed(i)
    X <- x0
    list.x <- c(X)
    list.y <- c()
    for(j in 1:T){
      X <- rho*X + rnorm(n=1, sd=sqrt(s2eps))
      Y <- X + rnorm(n=1, sd=sqrt(s2eta))
      list.x <- c(list.x, X)
      list.y <- c(list.y, Y)
    }
    
    dat <- matrix(data = list.y, nrow = 1)
    Z <- matrix(1, 1, 1)
    A <- matrix(0, 1, 1)
    R <- matrix(list("s2eta"), 1, 1)
    B <- matrix(list("rho"), 1, 1)
    U <- matrix(0, 1, 1)
    Q <- matrix(list("s2eps"), 1, 1)
    x0 <- matrix(x0, 1, 1)
    model.gen <- list(Z = Z, A = A, R = R, B = B, U = U,
                      Q = Q, x0 = x0, tinitx = 0)
    
    res <- MARSS(dat, 
                 model = model.gen, 
                 control=list(maxit=1000, conv.test.slope.tol=0.01),
                 silent = 1)
    list.rho <- c(list.rho, res$coef[2])
    list.s2eps <- c(list.s2eps, res$coef[3])
    list.s2eta <- c(list.s2eta, res$coef[1])
    Sys.sleep(0.5); setTxtProgressBar(pb, i/n)
  }
  return(data.frame(rho = list.rho, s2eps = list.s2eps, s2eta = list.s2eta))
}

simulationMAR1 <- function(n=1000, T=50, a11=-1.5, a12=0.7, a21=-4, a22=2, s2eps, s2eta, x01, x02){
  list.a11 <- c()
  list.a12 <- c()
  list.a21 <- c()
  list.a22 <- c()
  list.s2eps <- c()
  list.s2eta <- c()
  pb <- txtProgressBar(char=' \u27a4', style=3)
  for(i in 1:n){
    set.seed(i)
    X <- x0 <- c(x01,x02)
    list.y <- cbind()
    B <- rbind(c(a11, a12),
               c(a21, a22))
    for(t in 1:T){
      X <- B %*% X + rnorm(n=2, sd=sqrt(s2eps))
      Y <- X + rnorm(n=2, sd=sqrt(s2eta))
      list.y <- cbind(list.y, Y)
    }
    
    dat <- matrix(data = list.y, nrow = 2)
    Z <- diag(2)
    A <- matrix(0, 2, 1)
    R <- matrix(list(0), 2, 2)
    diag(R) <- c("s2eta", "s2eta")
    B <- matrix(list("a11", "a21", "a12", "a22"), 2, 2)
    U <- matrix(c(0, 0), 2, 1)
    Q <- matrix(list(0), 2, 2)
    diag(Q) <- c("s2eps", "s2eps")
    x0 <- matrix(c(x0), 2, 1)
    model.gen <- list(Z = Z, A = A, R = R, B = B, U = U,
                      Q = Q,
                      x0 = x0,
                      tinitx = 0)
    
    res <- MARSS(dat, 
                 model = model.gen, 
                 control=list(maxit=1000, conv.test.slope.tol=0.01),
                 silent = 1)
    list.a11 <- c(list.a11, res$coef[2])
    list.a12 <- c(list.a12, res$coef[4])
    list.a21 <- c(list.a21, res$coef[3])
    list.a22 <- c(list.a22, res$coef[5])
    list.s2eps <- c(list.s2eps, res$coef[6])
    list.s2eta <- c(list.s2eta, res$coef[1])
    
    Sys.sleep(0.5); setTxtProgressBar(pb, i/n)
  }
  
  return(data.frame(a11 = list.a11, a12 = list.a12,
                    a21 = list.a21, a22 = list.a22,
                    s2eps = list.s2eps, s2eta = list.s2eta))
  
}



## AR(1) with only stationary ----
# rho=0.45, s2eps=0.5, s2eta=0.5 

res1 <- simulationAR1(rho=0.45, s2eps=0.5, s2eta=0.5, x0=0)
save(res1, file = 'res1.RDATA')

round(mean(res1$rho), digits = 3) # 0.457
round(sd(res1$rho), digits = 3) # 0.291
round(mean(res1$s2eps), digits = 3) # 0.536
round(sd(res1$s2eps), digits = 3) # 0.369
round(mean(res1$s2eta), digits = 3) # 0.438
round(sd(res1$s2eta), digits = 3) # 0.330


# rho=0.45, s2eps=0.01, s2eta=0.1 

res2 <- simulationAR1(rho=0.45, s2eps=0.01, s2eta=0.1, x0=0)
save(res2, file = 'res2.RDATA')

round(mean(res2$rho), digits = 3) # 0.132
round(sd(res2$rho), digits = 3) # 0.425
round(mean(res2$s2eps), digits = 3) # 0.049
round(sd(res2$s2eps), digits = 3) # 0.034
round(mean(res2$s2eta), digits = 3) # 0.058
round(sd(res2$s2eta), digits = 3) # 0.034


# rho=0.7, s2eps=0.5, s2eta=0.5 

res3 <- simulationAR1(rho=0.7, s2eps=0.5, s2eta=0.5, x0=0)
save(res3, file = 'res3.RDATA')

round(mean(res3$rho), digits = 3) # 0.649
round(sd(res3$rho), digits = 3) # 0.205
round(mean(res3$s2eps), digits = 3) # 0.579
round(sd(res3$s2eps), digits = 3) # 0.348
round(mean(res3$s2eta), digits = 3) # 0.419
round(sd(res3$s2eta), digits = 3) # 0.292


# rho=0.7, s2eps=0.01, s2eta=0.1

res4 <- simulationAR1(rho=0.7, s2eps=0.01, s2eta=0.1, x0=0)
save(res4, file = 'res4.RDATA')

round(mean(res4$rho), digits = 3) # 0.298
round(sd(res4$rho), digits = 3) # 0.453
round(mean(res4$s2eps), digits = 3) # 0.044
round(sd(res4$s2eps), digits = 3) # 0.035
round(mean(res4$s2eta), digits = 3) # 0.065
round(sd(res4$s2eta), digits = 3) # 0.036


# rho=0, s2eps=0.5, s2eta=0.5 

res5 <- simulationAR1(rho=0, s2eps=0.5, s2eta=0.5, x0=0)
save(res5, file = 'res5.RDATA')

round(mean(res5$rho), digits = 3) # 0.063
round(sd(res5$rho), digits = 3) # 0.436
round(mean(res5$s2eps), digits = 3) # 0.412
round(sd(res5$s2eps), digits = 3) # 0.319
round(mean(res5$s2eta), digits = 3) # 0.543
round(sd(res5$s2eta), digits = 3) # 0.329


# rho=0, s2eps=0.01, s2eta=0.1 

res6 <- simulationAR1(rho=0, s2eps=0.01, s2eta=0.1, x0=0)
save(res6, file = 'res6.RDATA')

round(mean(res6$rho), digits = 3) # 0.025
round(sd(res6$rho), digits = 3) # 0.409
round(mean(res6$s2eps), digits = 3) # 0.051
round(sd(res6$s2eps), digits = 3) # 0.033
round(mean(res6$s2eta), digits = 3) # 0.055
round(sd(res6$s2eta), digits = 3) # 0.033


# Graphics ----
library(ggplot2)

load("/home/jgibaud/Documents/recherche/code/code ssm/res1.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res2.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res3.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res4.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res5.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res6.RDATA")

df <- rbind(res1[,2:3],
            res2[,2:3],
            res3[,2:3],
            res4[,2:3],
            res5[,2:3],
            res6[,2:3])

simu <- c(
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "rho *' = '* 0.7 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.7 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "rho *' = '* 0 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1"
)

df$simu <- as.factor(rep(simu, each=1000))

pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=3, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp

## AR(1) with transients: x0=5 ----
# rho=0.45, s2eps=0.5, s2eta=0.5 

res7 <- simulationAR1(rho=0.45, s2eps=0.5, s2eta=0.5, x0=5)
save(res7, file = 'res7.RDATA')

round(mean(res7$rho), digits = 3) # 0.449
round(sd(res7$rho), digits = 3) # 0.143
round(mean(res7$s2eps), digits = 3) # 0.485
round(sd(res7$s2eps), digits = 3) # 0.352
round(mean(res7$s2eta), digits = 3) # 0.492
round(sd(res7$s2eta), digits = 3) # 0.336


# rho=0.45, s2eps=0.01, s2eta=0.1 

res8 <- simulationAR1(rho=0.45, s2eps=0.01, s2eta=0.1, x0=5)
save(res8, file = 'res8.RDATA')

round(mean(res8$rho), digits = 3) # 0.448
round(sd(res8$rho), digits = 3) # 0.047
round(mean(res8$s2eps), digits = 3) # 0.016
round(sd(res8$s2eps), digits = 3) # 0.022
round(mean(res8$s2eta), digits = 3) # 0.091
round(sd(res8$s2eta), digits = 3) # 0.031


# rho=0.9, s2eps=0.5, s2eta=0.5

res9 <- simulationAR1(rho=0.9, s2eps=0.5, s2eta=0.5, x0=5)
save(res9, file = 'res9.RDATA')

round(mean(res9$rho), digits = 3) # 0.880
round(sd(res9$rho), digits = 3) # 0.059
round(mean(res9$s2eps), digits = 3) # 0.487
round(sd(res9$s2eps), digits = 3) # 0.266
round(mean(res9$s2eta), digits = 3) # 0.505
round(sd(res9$s2eta), digits = 3) # 0.238


# rho=0.9, s2eps=0.01, s2eta=0.1 

res10 <- simulationAR1(rho=0.9, s2eps=0.01, s2eta=0.1, x0=5)
save(res10, file = 'res10.RDATA')

round(mean(res10$rho), digits = 3) # 0.899
round(sd(res10$rho), digits = 3) # 0.010
round(mean(res10$s2eps), digits = 3) # 0.008
round(sd(res10$s2eps), digits = 3) # 0.008
round(mean(res10$s2eta), digits = 3) # 0.103
round(sd(res10$s2eta), digits = 3) # 0.025


## AR(1) with transients: x0=10 ----
# rho=0.45, s2eps=0.5, s2eta=0.5 

res11 <- simulationAR1(rho=0.45, s2eps=0.5, s2eta=0.5, x0=10)
save(res11, file = 'res11.RDATA')

round(mean(res11$rho), digits = 3) # 0.448
round(sd(res11$rho), digits = 3) # 0.080
round(mean(res11$s2eps), digits = 3) # 0.465
round(sd(res11$s2eps), digits = 3) # 0.322
round(mean(res11$s2eta), digits = 3) # 0.512
round(sd(res11$s2eta), digits = 3) # 0.320


# rho=0.45, s2eps=0.01, s2eta=0.1 

res12 <- simulationAR1(rho=0.45, s2eps=0.01, s2eta=0.1, x0=10)
save(res12, file = 'res12.RDATA')

round(mean(res12$rho), digits = 3) # 0.449
round(sd(res12$rho), digits = 3) # 0.023
round(mean(res12$s2eps), digits = 3) # 0.016
round(sd(res12$s2eps), digits = 3) # 0.021
round(mean(res12$s2eta), digits = 3) # 0.091
round(sd(res12$s2eta), digits = 3) # 0.031


# rho=0.9, s2eps=0.5, s2eta=0.5 

res13 <- simulationAR1(rho=0.9, s2eps=0.5, s2eta=0.5, x0=10)
save(res13, file = 'res13.RDATA')

round(mean(res13$rho), digits = 3) # 0.893
round(sd(res13$rho), digits = 3) # 0.031
round(mean(res13$s2eps), digits = 3) # 0.458
round(sd(res13$s2eps), digits = 3) # 0.244
round(mean(res13$s2eta), digits = 3) # 0.525
round(sd(res13$s2eta), digits = 3) # 0.232


# rho=0.9, s2eps=0.01, s2eta=0.1 

res14 <- simulationAR1(rho=0.9, s2eps=0.01, s2eta=0.1, x0=10)
save(res14, file = 'res14.RDATA')

round(mean(res14$rho), digits = 3) # 0.900
round(sd(res14$rho), digits = 3) # 0.005
round(mean(res14$s2eps), digits = 3) # 0.008
round(sd(res14$s2eps), digits = 3) # 0.008
round(mean(res14$s2eta), digits = 3) # 0.103
round(sd(res14$s2eta), digits = 3) # 0.025


# Graphics ----
library(ggplot2)

load("/home/jgibaud/Documents/recherche/code/code ssm/res7.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res8.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res9.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res10.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res11.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res12.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res13.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res14.RDATA")

df <- rbind(res7[,2:3],
            res8[,2:3],
            res9[,2:3],
            res10[,2:3])

simu <- c(
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "rho *' = '* 0.9 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.9 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1"
)


df$simu <- as.factor(rep(simu, each=1000))

pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=3, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp


df <- rbind(res11[,2:3],
            res12[,2:3],
            res13[,2:3],
            res14[,2:3])

simu <- c(
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.45 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "rho *' = '* 0.9 * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "rho *' = '* 0.9 * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1"
)

df$simu <- as.factor(rep(simu, each=1000))

pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=3, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp


# MAR(1) with only stationary ----
# s2eps=0.5, s2eta=0.5 

res15 <- simulationMAR1(s2eps=0.5, s2eta=0.5, x01=0, x02=0)
save(res15, file = 'res15.RDATA')

round(mean(res15$a11), digits = 3) # -1.541
round(sd(res15$a11), digits = 3) # 0.313
round(mean(res15$a12), digits = 3) # 0.724
round(sd(res15$a12), digits = 3) # 0.143
round(mean(res15$a21), digits = 3) # -4.097
round(sd(res15$a21), digits = 3) # 0.786
round(mean(res15$a22), digits = 3) # 2.027
round(sd(res15$a22), digits = 3) # 0.288
round(mean(res15$s2eps), digits = 3) # 0.495
round(sd(res15$s2eps), digits = 3) # 0.180
round(mean(res15$s2eta), digits = 3) # 0.475
round(sd(res15$s2eta), digits = 3) # 0.108

# s2eps=0.01, s2eta=0.1 

res16 <- simulationMAR1(s2eps=0.01, s2eta=0.1, x01=0, x02=0)
save(res16, file = 'res16.RDATA')

round(mean(res16$a11), digits = 3) # -1.090
round(sd(res16$a11), digits = 3) # 1.055
round(mean(res16$a12), digits = 3) # 0.599
round(sd(res16$a12), digits = 3) # 0.450
round(mean(res16$a21), digits = 3) # -2.892
round(sd(res16$a21), digits = 3) # 2.908
round(mean(res16$a22), digits = 3) # 1.640
round(sd(res16$a22), digits = 3) # 0.980
round(mean(res16$s2eps), digits = 3) # 0.017
round(sd(res16$s2eps), digits = 3) # 0.02
round(mean(res16$s2eta), digits = 3) # 0.094
round(sd(res16$s2eta), digits = 3) # 0.023

# s2eps=0.01, s2eta=1 

res21 <- simulationMAR1(s2eps=0.01, s2eta=1, x01=0, x02=0)
save(res21, file = 'res21.RDATA')

round(mean(res21$a11), digits = 3) # -0.006
round(sd(res21$a11), digits = 3) # 0.789
round(mean(res21$a12), digits = 3) # 0.177
round(sd(res21$a12), digits = 3) # 0.750
round(mean(res21$a21), digits = 3) # -0.100
round(sd(res21$a21), digits = 3) # 1.709
round(mean(res21$a22), digits = 3) # 0.548
round(sd(res21$a22), digits = 3) # 0.741
round(mean(res21$s2eps), digits = 3) # 0.217
round(sd(res21$s2eps), digits = 3) # 0.306
round(mean(res21$s2eta), digits = 3) # 0.739
round(sd(res21$s2eta), digits = 3) # 0.287

# MAR(1) with x0=(-5,5) ----
# s2eps=0.5, s2eta=0.5 

res17 <- simulationMAR1(s2eps=0.5, s2eta=0.5, x01=-5, x02=5)
save(res17, file = 'res17.RDATA')

round(mean(res17$a11), digits = 3) # -1.504
round(sd(res17$a11), digits = 3) # 0.102
round(mean(res17$a12), digits = 3) # 0.706
round(sd(res17$a12), digits = 3) # 0.056
round(mean(res17$a21), digits = 3) # -4.000
round(sd(res17$a21), digits = 3) # 0.162
round(mean(res17$a22), digits = 3) # 2.001
round(sd(res17$a22), digits = 3) # 0.061
round(mean(res17$s2eps), digits = 3) # 0.477
round(sd(res17$s2eps), digits = 3) # 0.110
round(mean(res17$s2eta), digits = 3) # 0.483
round(sd(res17$s2eta), digits = 3) # 0.107

# s2eps=0.01, s2eta=0.1

res18 <- simulationMAR1(s2eps=0.01, s2eta=0.1, x01=-5, x02=5)
save(res18, file = 'res18.RDATA')

round(mean(res18$a11), digits = 3) # -1.502
round(sd(res18$a11), digits = 3) # 0.031
round(mean(res18$a12), digits = 3) # 0.703
round(sd(res18$a12), digits = 3) # 0.025
round(mean(res18$a21), digits = 3) # -3.996
round(sd(res18$a21), digits = 3) # 0.063
round(mean(res18$a22), digits = 3) # 2.002
round(sd(res18$a22), digits = 3) # 0.025
round(mean(res18$s2eps), digits = 3) # 0.009
round(sd(res18$s2eps), digits = 3) # 0.003
round(mean(res18$s2eta), digits = 3) # 0.097
round(sd(res18$s2eta), digits = 3) # 0.021

# s2eps=0.01, s2eta=1 

res22 <- simulationMAR1(s2eps=0.01, s2eta=1, x01=-5, x02=5)
save(res22, file = 'res22.RDATA')

round(mean(res22$a11), digits = 3) # -1.530
round(sd(res22$a11), digits = 3) # 0.090
round(mean(res22$a12), digits = 3) # 0.755
round(sd(res22$a12), digits = 3) # 0.080
round(mean(res22$a21), digits = 3) # -3.885
round(sd(res22$a21), digits = 3) # 0.198
round(mean(res22$a22), digits = 3) # 2.037
round(sd(res22$a22), digits = 3) # 0.075
round(mean(res22$s2eps), digits = 3) # 0.009
round(sd(res22$s2eps), digits = 3) # 0.008
round(mean(res22$s2eta), digits = 3) # 0.974
round(sd(res22$s2eta), digits = 3) # 0.178


# MAR(1) with x0=(-10,10) ----
# s2eps=0.5, s2eta=0.5 

res19 <- simulationMAR1(s2eps=0.5, s2eta=0.5, x01=-10, x02=10)
save(res19, file = 'res19.RDATA')

round(mean(res19$a11), digits = 3) # -1.502
round(sd(res19$a11), digits = 3) # 0.058
round(mean(res19$a12), digits = 3) # 0.703
round(sd(res19$a12), digits = 3) # 0.032
round(mean(res19$a21), digits = 3) # -3.999
round(sd(res19$a21), digits = 3) # 0.086
round(mean(res19$a22), digits = 3) # 2.001
round(sd(res19$a22), digits = 3) # 0.032
round(mean(res19$s2eps), digits = 3) # 0.476
round(sd(res19$s2eps), digits = 3) # 0.106
round(mean(res19$s2eta), digits = 3) # 0.484
round(sd(res19$s2eta), digits = 3) # 0.108

# s2eps=0.01, s2eta=0.1

res20 <- simulationMAR1(s2eps=0.01, s2eta=0.1, x01=-10, x02=10)
save(res20, file = 'res20.RDATA')

round(mean(res20$a11), digits = 3) # -1.501
round(sd(res20$a11), digits = 3) # 0.016
round(mean(res20$a12), digits = 3) # 0.702
round(sd(res20$a12), digits = 3) # 0.013
round(mean(res20$a21), digits = 3) # -3.997
round(sd(res20$a21), digits = 3) # 0.032
round(mean(res20$a22), digits = 3) # 2.001
round(sd(res20$a22), digits = 3) # 0.012
round(mean(res20$s2eps), digits = 3) # 0.009
round(sd(res20$s2eps), digits = 3) # 0.003
round(mean(res20$s2eta), digits = 3) # 0.097
round(sd(res20$s2eta), digits = 3) # 0.021

# s2eps=0.01, s2eta=1

res23 <- simulationMAR1(s2eps=0.01, s2eta=1, x01=-10, x02=10)
save(res23, file = 'res23.RDATA')

round(mean(res23$a11), digits = 3) # -1.518
round(sd(res23$a11), digits = 3) # 0.045
round(mean(res23$a12), digits = 3) # 0.729
round(sd(res23$a12), digits = 3) # 0.040
round(mean(res23$a21), digits = 3) # -3.945
round(sd(res23$a21), digits = 3) # 0.099
round(mean(res23$a22), digits = 3) # 2.021
round(sd(res23$a22), digits = 3) # 0.038
round(mean(res23$s2eps), digits = 3) # 0.009
round(sd(res23$s2eps), digits = 3) # 0.008
round(mean(res23$s2eta), digits = 3) # 0.973
round(sd(res23$s2eta), digits = 3) # 0.177

#Graphics ----
library(ggplot2)

load("/home/jgibaud/Documents/recherche/code/code ssm/res15.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res16.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res17.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res18.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res19.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res20.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res21.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res22.RDATA")
load("/home/jgibaud/Documents/recherche/code/code ssm/res23.RDATA")


df <- rbind(res15[,5:6],
            res16[,5:6],
            res21[,5:6])


simu <- c(
  "x[0] *' = '* (0*','*0) * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "x[0] *' = '* (0*','*0) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "x[0] *' = '* (0*','*0) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 1"
)

df$simu <- as.factor(rep(simu, each=1000))


# pdf(file = "s2etaWRTs2epsMAR(1).pdf", width = 7, height = 10)
pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=2, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp
# dev.off()


df <- rbind(res17[,5:6],
            res18[,5:6],
            res22[,5:6])


simu <- c(
  "x[0] *' = '* (-5*','*5) * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "x[0] *' = '* (-5*','*5) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "x[0] *' = '* (-5*','*5) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 1"
)

df$simu <- as.factor(rep(simu, each=1000))


# pdf(file = "s2etaWRTs2epsMAR(1).pdf", width = 7, height = 10)
pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=5, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp
# dev.off()


df <- rbind(res19[,5:6],
            res20[,5:6],
            res23[,5:6])


simu <- c(
  "x[0] *' = '* (-10*','*10) * ', ' * sigma[epsilon]^2 *' = '* 0.5 * ', ' * sigma[eta]^2 *' = ' * 0.5",
  "x[0] *' = '* (-10*','*10) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 0.1",
  "x[0] *' = '* (-10*','*10) * ', ' * sigma[epsilon]^2 *' = '* 0.01 * ', ' * sigma[eta]^2 *' = ' * 1"
)

df$simu <- as.factor(rep(simu, each=1000))


# pdf(file = "s2etaWRTs2epsMAR(1).pdf", width = 7, height = 10)
pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=5, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp
# dev.off()



# Case study ----

#define functions for parametric bootstrap
BootAR1 <- function(n=1000, T=35, rho, s2eps, s2eta, x0){
  list.x0 <- c()
  list.rho <- c()
  list.s2eps <- c()
  list.s2eta <- c()
  pb <- txtProgressBar(char=' \u27a4', style=3)
  for(i in 1:n){
    set.seed(i)
    X <- x0
    list.x <- c(X)
    list.y <- c()
    for(j in 1:T){
      X <- rho*X + rnorm(n=1, sd=sqrt(s2eps))
      Y <- X + rnorm(n=1, sd=sqrt(s2eta))
      list.x <- c(list.x, X)
      list.y <- c(list.y, Y)
    }
    
    dat <- matrix(data = list.y, nrow = 1)
    Z <- matrix(1, 1, 1)
    A <- matrix(0, 1, 1)
    R <- matrix(list("s2eta"), 1, 1)
    B <- matrix(list("rho"), 1, 1)
    U <- matrix(0, 1, 1)
    Q <- matrix(list("s2eps"), 1, 1)
    model.gen <- list(Z = Z, A = A, R = R, B = B, U = U,
                      Q = Q, tinitx = 0)
    
    res <- MARSS(dat, 
                 model = model.gen, 
                 control=list(maxit=1000, conv.test.slope.tol=0.01),
                 silent = 1)
    list.rho <- c(list.rho, res$coef[2])
    list.s2eps <- c(list.s2eps, res$coef[3])
    list.s2eta <- c(list.s2eta, res$coef[1])
    list.x0 <- c(list.x0, res$coef[4])
    Sys.sleep(0.5); setTxtProgressBar(pb, i/n)
  }
  return(data.frame(rho = list.rho, s2eps = list.s2eps, s2eta = list.s2eta, x0 = list.x0))
}

BootMAR1 <- function(n=1000, T=35, a11, a12, a21, a22, s2eps, s2eta, x01, x02){
  list.a11 <- c()
  list.a12 <- c()
  list.a21 <- c()
  list.a22 <- c()
  list.s2eps <- c()
  list.s2eta <- c()
  list.x01 <- c()
  list.x02 <- c()
  pb <- txtProgressBar(char=' \u27a4', style=3)
  for(i in 1:n){
    set.seed(i)
    X <- x0 <- c(x01,x02)
    list.y <- cbind()
    B <- rbind(c(a11, a12),
               c(a21, a22))
    for(t in 1:T){
      X <- B %*% X + rnorm(n=2, sd=sqrt(s2eps))
      Y <- X + rnorm(n=2, sd=sqrt(s2eta))
      list.y <- cbind(list.y, Y)
    }
    
    dat <- matrix(data = list.y, nrow = 2)
    Z <- diag(2)
    A <- matrix(0, 2, 1)
    R <- matrix(list(0), 2, 2)
    diag(R) <- c("s2eta", "s2eta")
    B <- matrix(list("a11", "a21", "a12", "a22"), 2, 2)
    U <- matrix(c(0, 0), 2, 1)
    Q <- matrix(list(0), 2, 2)
    diag(Q) <- c("s2eps", "s2eps")
    model.gen <- list(Z = Z, A = A, R = R, B = B, U = U,
                      Q = Q, tinitx = 0)
    
    res <- MARSS(dat, 
                 model = model.gen, 
                 control=list(maxit=1000, conv.test.slope.tol=0.01),
                 silent = 1)
    list.a11 <- c(list.a11, res$coef[2])
    list.a12 <- c(list.a12, res$coef[4])
    list.a21 <- c(list.a21, res$coef[3])
    list.a22 <- c(list.a22, res$coef[5])
    list.s2eps <- c(list.s2eps, res$coef[6])
    list.s2eta <- c(list.s2eta, res$coef[1])
    list.x01 <- c(list.x01, res$coef[7])
    list.x02 <- c(list.x02, res$coef[8])
    
    Sys.sleep(0.5); setTxtProgressBar(pb, i/n)
  }
  
  return(data.frame(a11 = list.a11, a12 = list.a12,
                    a21 = list.a21, a22 = list.a22,
                    s2eps = list.s2eps, s2eta = list.s2eta,
                    x01 = list.x01, x02 = list.x02))
  
}

#load data
data <- read.csv2("/home/jgibaud/Documents/recherche/code/code ssm/jzo12716-sup-0002-datas2.csv", header = T)

# Graphics
df <- as.data.frame(cbind(rep(data$year_N, times = 2),
                          c(data$N_deer_adj_SNP, data$N_chamois_adj),
                          rep(c("Deer", "Chamois"), each = length(data$year_N))))
colnames(df) <- c("Year", "Abundance", "Taxa")
df$Taxa <- as.factor(df$Taxa)
df$Year <- as.numeric(df$Year)
df$Abundance <- as.numeric(df$Abundance)


pp <- ggplot(data=df, aes(x=Year, y=Abundance, colour = Taxa)) + geom_line() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25),
        legend.text=element_text(size = 25), legend.title=element_text(size = 20))+
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015))

# pdf(file = "/home/jgibaud/Documents/recherche/code/code ssm/Abundance.pdf", width = 12, height = 7.5)
pp
# dev.off()


df$Abundance <- log(df$Abundance)
pp <- ggplot(data=df, aes(x=Year, y=Abundance, colour = Taxa)) + geom_line()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25),
        legend.text=element_text(size = 18), legend.title=element_text(size = 15))+
  scale_x_continuous(breaks=c(1985, 1990, 1995, 2000, 2005, 2010, 2015))

# pdf(file = "/home/jgibaud/Documents/recherche/code/code ssm/LogAbundance.pdf", width = 12, height = 7.5)
pp
# dev.off()

#run MARSS for deer
Z <- matrix(1, 1, 1)
A <- matrix(0, 1, 1)
R <- matrix(list("s2eta"), 1, 1)
B <- matrix(list("rho"), 1, 1)
U <- matrix(0, 1, 1)
Q <- matrix(list("s2eps"), 1, 1)

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, tinitx = 0)
control <- list(maxit=50000, conv.test.slope.tol=0.01, abstol=0.01)

Y.deer <- log(as.numeric(data$N_deer_adj_SNP))
res.deer <- MARSS(Y.deer,
                  model = model.gen,
                  control = control)
res.deer

boot.deer <- BootAR1(rho = res.deer$par$B,
                     s2eps = res.deer$par$Q,
                     s2eta = res.deer$par$R,
                     x0 = res.deer$par$x0)
boot.deer


#run MARSS for chamois
Y.chamois <- log(as.numeric(data$N_chamois_adj))
res.chamois <- MARSS(Y.chamois,
                     model = model.gen,
                     control = control)
res.chamois

boot.chamois <- BootAR1(rho = res.chamois$par$B,
                        s2eps = res.chamois$par$Q,
                        s2eta = res.chamois$par$R,
                        x0 = res.chamois$par$x0)
boot.chamois


#run MARSS for ungulate
Y.ungulate <- matrix(data = rbind(Y.deer, Y.chamois), nrow = 2)

Z <- diag(2)
A <- matrix(0, 2, 1)
R <- matrix(list(0), 2, 2)
diag(R) <- c("s2eta", "s2eta")
B <- matrix(list("a11", "a21", "a12", "a22"), 2, 2)
U <- matrix(c(0, 0), 2, 1)
Q <- matrix(list(0), 2, 2)
diag(Q) <- c("s2eps", "s2eps")
model.ungulate <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, tinitx = 0)

res.ungulate <- MARSS(Y.ungulate,
                      model = model.ungulate,
                      control = control)
res.ungulate

boot.ungulate <- BootMAR1(a11 = res.ungulate$coef[2],
                          a12 = res.ungulate$coef[4], 
                          a21 = res.ungulate$coef[3],
                          a22 = res.ungulate$coef[5],
                          s2eps = res.ungulate$coef[6],
                          s2eta = res.ungulate$coef[1],
                          x01 = res.ungulate$coef[7],
                          x02 = res.ungulate$coef[8])
boot.ungulate

save(res.chamois,
     res.deer,
     res.ungulate,
     boot.chamois,
     boot.deer,
     boot.ungulate,
     file = "/home/jgibaud/Documents/recherche/code/code ssm/CaseStudy.RDATA")

#Graphics----

load("/home/jgibaud/Documents/recherche/code/code ssm/CaseStudy.RDATA")

df <- rbind(boot.chamois[,2:3],
            boot.deer[,2:3],
            boot.ungulate[,5:6])


df$simu <- as.factor(c(rep("Chamois", 1000),
                       rep("Deer", 1000),
                       rep("Ungulates", 998)))


pdf(file = "/home/jgibaud/Documents/recherche/code/code ssm/VariancesCaseStudy.pdf", width = 7, height = 10)
pp <- ggplot(df, aes(s2eps, s2eta)) + geom_point(size=0.5)
pp <- pp + facet_wrap(~ simu, scales = "free", nrow=2, ncol=2, labeller = label_parsed)
pp <- pp + xlab("State variance") + ylab("Observation variance")
pp
dev.off()
