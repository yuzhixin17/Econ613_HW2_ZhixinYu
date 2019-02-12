############################# R Script for Econ613_HW#2 #############################
############################# Feb 9,2019 #############################


######################## Exercise 1  Data Creation ########################

# setting a seed
set.seed(10)

# create objects X1,X2,X3 and eps
# generating random numbers from specific distributions
X1 <- runif(n = 10000, min = 1, max = 3)
X2 <- rgamma(n = 10000, shape = 3, rate = 0.5)
X3 <- rbinom(n = 10000, size = 1, prob = 0.3)
eps <- rnorm(n = 10000, mean = 2, sd = 1)

# create variables Y
Y <- 0.5+1.2*X1+(-0.9*X2)+0.1*X3+eps

# create dummy variable ydum
ydum <- as.numeric(Y > mean(Y))


# save created data as csv.file
datset <- cbind(1,X1,X2,X3,Y,ydum)
write.csv(datset,"E1_Data creation")

######################## Exercise 2  OLS ########################

## Question 1 ##
## Calulate the correlation between Y and X1
Cor_XY <- (10000*(sum(X1*Y))-(sum(X1)*sum(Y)))/sqrt((10000*sum(X1^2)-(sum(X1))^2)*(10000*sum(Y^2)-(sum(Y))^2))
Cor_XY
# [1] 0.2249412  ##it is quite different from 1.2 

## Question 2 ##
## We're interested in the outcome of the regression of Y on X where X=(1,X1,X2,X3)
X <- as.matrix(cbind(1,X1,X2,X3))
y <- as.matrix(Y)

## Question 3 ##
## Estimate the coefficients beta using formula : beta_hat = ((X'X)^(-1))X'y
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%y
beta_hat
# [,1]
#    2.5120783
#X1  1.1749022
#X2 -0.8967416
#X3  0.1091988


## Question 4 ##
## Calculate the standard errors

# (1) Using the standard formulas of the OLS
# residuals
resid <- y - X%*%beta_hat
# calculate the OLS estimate for sigma square
sigma2_hat <- (t(resid)%*%resid)/(nrow(X)-ncol(X))
# estimate of V[beta_hat|X], which is conditional covariance matrix of the least squares slope estimator
vcov_beta_hat <- c(sigma2_hat)*solve(t(X)%*%X)
# estimate of standard errors
sqrt(diag(vcov_beta_hat)) 
#X1          X2          X3 
#0.040312946 0.017053768 0.002884871 0.021565945 

# (2) Using bootstrap with 49 replications
# create a dataset combining X and Y
datset1 <- cbind(X,Y)

# create an empty matrix boot49
boot49 <- NULL
# write a for loop to calculate each of 49 sample's standard errors
for(i in 1:49){
    # sample rows from datset with replecament
    dat49 <- datset1[sample(nrow(datset1), 10000, replace = TRUE),]
    # form new X
    X49 <- as.matrix(cbind(1,dat49[,c(2:4)]))
    # form new Y
    y49 <- as.matrix(dat49[,5])
    # calculate se for each of 49 samples
    beta_hat49 <- solve(t(X49)%*%X49)%*%t(X49)%*%y49
    resid49 <- y49 - X49%*%beta_hat49
    sigma2_hat49 <- (t(resid49)%*%resid49)/(nrow(X49)-ncol(X49))
    vcov_beta_hat49 <- c(sigma2_hat49)*solve(t(X49)%*%X49)
    se49 <- sqrt(diag(vcov_beta_hat49))
    # save results in the matrix boot49
    boot49 <- rbind(boot49,se49)
    }

# Calculate the average of standard error
colMeans(boot49)
#X1          X2          X3 
#0.040310899 0.017059116 0.002891144 0.021569302 

# (2) Using bootstrap with 499 replications

# create an empty matrix boot499
boot499 <- NULL
# write a for loop to calculate each of 499 sample's standard errors
for(i in 1:499){
    # sample rows from datset with replecament
    dat499 <- datset1[sample(nrow(datset1), 10000, replace = TRUE),]
    # form new X
    X499 <- as.matrix(cbind(1,dat499[,c(2:4)]))
    # form new Y
    y499 <- as.matrix(dat499[,5])
    # calculate se for each of 49 samples
    beta_hat499 <- solve(t(X499)%*%X499)%*%t(X499)%*%y499
    resid499 <- y499 - X499%*%beta_hat499
    sigma2_hat499 <- (t(resid499)%*%resid499)/(nrow(X499)-ncol(X499))
    vcov_beta_hat499 <- c(sigma2_hat499)*solve(t(X499)%*%X499)
    se499 <- sqrt(diag(vcov_beta_hat499))
    # save results in the matrix boot49
    boot499 <- rbind(boot499,se499)
    }

# Calculate the average of standard error
colMeans(boot499)
#X1          X2          X3 
#0.040308007 0.017056341 0.002884332 0.021562551



######################## Exercise 3  Numerical Optinization ########################

## Question 1 ##
# write a function that returns the likelihood of the probit
Logl_probit <- function(beta, ydum, X){
    logl = sum(ydum*log(pnorm(X%*%beta))+(1-ydum)*log(1-pnorm(X%*%beta)))
    logl = -logl
    return(logl)
}

## Question 2 ##
# implement the steepest ascent optimization algorithm to maximize that likelihood
# write a function calculating gradients
grad <- function(beta, ydum, X){
    n = length(ydum)           # sample size
    k = length(beta)           # number of coefficients
    g = t(matrix(rep(dnorm(X%*%beta)/(pnorm(X%*%beta)),k),nrow=n)*X) %*% ydum - 
        t(matrix(rep(dnorm(X%*%beta)/(1-pnorm(X%*%beta)),k),nrow=n)*X) %*% (1-ydum)
    g = -g
    return(g)
}
# set starting value and alpha value
beta_old = c(-1,-1,-1,-1)
beta_new = c(0,0,0,0)
alpha = 0.00001
# implement the steepest ascent optimization algorithm and set tolerance = 0.00000001
while (abs(beta_new - beta_old) > 0.00000001){
      beta_old <- beta_new
      beta_new <- beta_old - alpha*grad(beta_old, ydum, X)    
}
print(beta_new)
#[,1]
#    3.0570853
#X1  1.1857779
#X2 -0.9113114
#X3  0.1075399

## Question 3 ##
# How different are the parameter from the true parameter
# use probit package to report true parameter
myprobit <- glm(ydum~1+X1+X2+X3,family = binomial(link = "probit"))
T_beta <- summary(myprobit)$coefficients[,1]
T_beta
#(Intercept)          X1          X2          X3 
#  3.0570959   1.1857762  -0.9113125   0.1075392

# save true parameter and our numerical optimization parameter
E3_parameter <- cbind(T_beta,beta_new)
colnames(E3_parameter) <- c("T_beta","E_beta")
# we can find that the difference between the true parameter and our optimized parameter is very small
# save result
write.csv(E3_parameter,"E3")
 


######################## Exercise 4  Disecrete Choice ########################

# write and optimize probit, logit and linear probability model 
# write log-likelihood function for probit model
Probit <- function(beta, ydum, X){
    phi = pnorm(X%*%beta)
    logl = sum(ydum*log(phi)+(1-ydum)*log(1-phi))
    logl = -logl
    return(logl)
}
# write log-likelihood function for logit model
Logit <- function(beta, ydum, X){
    gam = exp(X%*%beta)/(1+exp(X%*%beta))
    logl = sum(ydum*log(gam)+(1-ydum)*log(1-gam))
    logl = -logl
    return(logl)
}
# calculating mean square for linear model
Linear <- function(beta, ydum, X){
    min<-t(ydum-X%*%beta)%*%(ydum-X%*%beta)
    return(min)
}

# optimization
optim(c(0,0,0,0), Probit, ydum = ydum, X=X)$par
#[1]  3.0567133  1.1856985 -0.9112391  0.1077866
optim(c(0,0,0,0), Logit, ydum = ydum, X=X)$par
#[1]  5.4608550  2.1159251 -1.6269076  0.1763368
optim(c(0,0,0,0), Linear, ydum = ydum, X=X)$par
#[1]  0.90031483  0.14240235 -0.10499390  0.00937559


# comparison
# we may find out that the difference between these two methods are very small
summary(myprobit)$coefficients[,1]
#(Intercept)          X1          X2          X3 
#  3.0570959   1.1857762  -0.9113125   0.1075392 
summary(glm(ydum~1+X1+X2+X3,family = binomial(link = "logit")))$coefficients[,1]
#(Intercept)          X1          X2          X3 
#  5.4597311   2.1159724  -1.6267401   0.1761663  
summary(lm(ydum~1+X1+X2+X3))$coefficients[,1]
#(Intercept)           X1           X2           X3 
#  0.899926742  0.142498386 -0.104961892  0.009316183 

# interpretation
# we can interprete that X1 and X3 have positive effect, while X2 has negative effect on Y.
# they are all significant, except for X3 in linear model
# for those are significant, they are siginificant at ‘***’ 0.001 level, except for vector X3 in probit and logit model are siginificant at ‘*’ 0.05 level
summary(myprobit)
#Call:
#glm(formula = ydum ~ 1 + X1 + X2 + X3, family = binomial(link = "probit"))

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-3.6510  -0.1085   0.0068   0.2404   3.2485  

#Coefficients:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  3.05710    0.10076  30.340   <2e-16 ***
#X1           1.18578    0.04355  27.228   <2e-16 ***
#X2          -0.91131    0.01864 -48.896   <2e-16 ***
#X3           0.10754    0.04760   2.259   0.0239 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(glm(ydum~1+X1+X2+X3,family = binomial(link = "logit")))
#Call:
#glm(formula = ydum ~ 1 + X1 + X2 + X3, family = binomial(link = "logit"))

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-3.2838  -0.1489   0.0374   0.2568   3.0310  

#Coefficients:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  5.45973    0.18840  28.979   <2e-16 ***
#X1           2.11597    0.08063  26.242   <2e-16 ***
#X2          -1.62674    0.03687 -44.123   <2e-16 ***
#X3           0.17617    0.08533   2.065    0.039 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(lm(ydum~1+X1+X2+X3))
#Call:
#lm(formula = ydum ~ 1 + X1 + X2 + X3)

#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.91032 -0.26830  0.05693  0.24758  1.64454 

#Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)  0.8999267  0.0134312   67.003   <2e-16 ***
#X1           0.1424984  0.0056819   25.079   <2e-16 ***
#X2          -0.1049619  0.0009612 -109.203   <2e-16 ***
#X3           0.0093162  0.0071852    1.297    0.195    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




######################## Exercise 5  Marginal Effect ########################

## Question 1 ##
# (1) compute the marginal effects of X and Y according to probit model
# Here we calculate Mean level of Marginal Effect

# (a) Probit Model

# run probit model
myprobit <- glm(ydum~1+X1+X2+X3,family = binomial(link = "probit"))
# save betas from probit model
beta_probit <- as.matrix(coef(myprobit))
# calculate mean of X
X_mean <- as.matrix(colMeans(X))
# transpose beta_probit for later calculation
beta_probit <- t(beta_probit)
# write a function calculating marginal effect of probit model
ME_probit <- function(beta){
    f = dnorm(X_mean%*%beta)%*%t(beta)
    return(f)
}

# (b) logit Model

# run logit model
mylogit <- glm(ydum~1+X1+X2+X3,family = binomial(link = "logit"))
# save betas from logit model
beta_logit <- as.matrix(coef(mylogit))
# transpose beta_logit for later calculation
beta_logit <- t(beta_logit)
# write a function calculating marginal effect of logit model
ME_logit <- function(beta){
    f = exp(X_mean%*%beta)%*%t(beta)/(1+(exp(X_mean%*%beta)%*%t(beta)))^2
    return(f)
}

## Question 2 ##

# (1) compute the standard deviations using the delta method
# Probit Model
# load package
install.packages("C:/Users/Administrator/Downloads/numDeriv_2016.8-1.zip", repos = NULL, type = "win.binary")
library(numDeriv)
# calculate jacobian matrix of marginal effect
jac_p <- as.matrix(jacobian(ME_probit, beta_probit))
# extract variance-covariance matrix
vcov_probit <- vcov(myprobit)
# use delta method do the calculation
delta_probit <- t(jac_p)%*%vcov_probit%*%jac_p
# calculate standard deviation
sqrt(diag(delta_probit))
#[1] 0.003821026 0.018132841 0.019075044 0.037518302

# Logit Model
jac_l <- as.matrix(jacobian(ME_logit, beta_logit))
vcov_logit <- vcov(mylogit)
delta_logit <- t(jac_l)%*%vcov_logit%*%jac_l
sqrt(diag(delta_logit))
#[1] 1.036827e-03 2.347287e-04 2.367264e-05 8.419558e-05


# (2) compute the standard deviations using bootstrap

# (a) Probit Model
boot_p499 <- NULL
# write a for loop to calculate each of 499 sample's standard deviation
for(i in 1:499){
    # sample rows from datset with replecament
    dat_p_499 <- as.data.frame(datset1[sample(nrow(datset1), 10000, replace = TRUE),])
    # form new X
    X_p499 <- as.matrix(dat_p_499[,c(1:4)])
    # calculate new mean of X
    X_mean_p499 <- as.matrix(colMeans(X_p499))
    # rewrite new marginal effect function
    ME_probit499 <- function(beta){
    f = dnorm(X_mean_p499%*%beta)%*%t(beta)
    return(f)
    }
    # calculate standard deviation for each of 499 sample
    jac_p499 <- as.matrix(jacobian(ME_probit499, beta_probit))
    vcov_probit499 <- vcov(myprobit)
    delta_probit499 <- t(jac_p)%*%vcov_probit499%*%jac_p
    sd499 <- sqrt(diag(delta_probit499))
    # save results in the matrix boot49
    boot_p499 <- rbind(boot_p499,sd499)
    }
colMeans(boot_p499) 
#[1] 0.003821026 0.018132841 0.019075044 0.037518302


# (b) Logit Model
boot_l499 <- NULL
# write a for loop to calculate each of 499 sample's standard deviation
for(i in 1:499){
    # sample rows from datset with replecament
    dat_l_499 <- as.data.frame(datset1[sample(nrow(datset1), 10000, replace = TRUE),])
    # form new X
    X_l499 <- as.matrix(dat_l_499[,c(1:4)])
    # calculate new mean of X
    X_mean_l499 <- as.matrix(colMeans(X_l499))
    # rewrite new marginal effect function
    ME_logit499 <- function(beta){
    f = exp(X_mean_l499%*%beta)%*%t(beta)/(1+(exp(X_mean_l499%*%beta)%*%t(beta)))^2
    return(f)
    }
    # calculate standard deviation for each of 499 sample
    jac_l499 <- as.matrix(jacobian(ME_logit499, beta_logit))
    vcov_logit499 <- vcov(mylogit)
    delta_logit499 <- t(jac_l)%*%vcov_logit499%*%jac_l
    sd499 <- sqrt(diag(delta_logit499))
    # save results in the matrix boot49
    boot_l499 <- rbind(boot_l499,sd499)
    }
colMeans(boot_l499) 
#[1] 1.036827e-03 2.347287e-04 2.367264e-05 8.419558e-05





########## Thank you for reading! Feel free making comments! Have a good day! ##########


##################################### The End ###############################################

