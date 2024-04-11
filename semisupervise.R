# Load packages

library(mev)
library(ggplot2)
library(sf)
library(rgdal)
library(rnaturalearth)
library(dplyr)
library(broom)

# Convert the first n observations by using the empirical
# distribution function based on all m observations 
tilde <- function (y, n){
    m <- length(y)
    y1 <- -log(1-rank(y)/ m+(1/(2*m)))
    return(y1[1:n])
}

# Improved gamma, with g=0
# data should be a two-column dataframe
# the first column, first n row, provides the main variable to be considered
# the second column, provides the covariate with more observations
# the first n cases match as pairs
# n is the sample size for the main variable
# k is a single value indicating the value k to be used

imp_gamma_mle <- function(data,n,k, k1){
    # Obtain m
    m <- dim(data)[1] - n
    
    # Obtain x
    x <- data[1:n,1]
    x_rank <- rank(x)

    # Prepare y_tilde
    y_tilde <- tilde(data[,2],n)
    y_rank <- rank(y_tilde)
    y_origin <-data[,2]

    # Get the thresholds
    x_threshold=sort(x)[n-k]
    y_threshold=sort(y_tilde)[n-k]
    y_origin_threshold=sort(y_origin)[m+n-k1]

    # MLE
    ml <- gp.fit(x,x_threshold,method="zhang") 
    gamma_1 <- ml$approx.mean[2] 
    sigma_1 <- ml$approx.mean[1] 
    var_gamma <- (1+gamma_1)^2
    
    ml1 <- gp.fit(y_tilde,y_threshold,method="zhang") 
    gamma_2 <- ml1$approx.mean[2]
    sigma_2 <- ml1$approx.mean[1]

    ml2 <- gp.fit(y_origin,y_origin_threshold,method="zhang") 
    gamma_origin <- ml2$approx.mean[2] 
    sigma_origin <- ml2$approx.mean[1] 
    var_gamma_origin <- (1+gamma_origin)^2

    # Calculate R11 and two integrals
    select <- (x>x_threshold & y_tilde>y_threshold)
    # R11 is the estimator of R(1,1)
    R11 <- sum(select)/k

    # Rint1 is the estimator of the integral from 0 to 1 on R(s,1)/s^{1-gamma}
    # The rationale is that R(s,1) is estimated by counting 
    # x_rank >n- [ks] AND y_rank >n-k
    # This only occurs fo observations with in "select" (x_rank > n-k)
    # For those observations, s > (n-x_rank+1)/k
    # The integral of s^{gamma-1} is then s^{gamma}/gamma
    # So we get the integral as 1-((n-x_rank+1)/k)^gamma
    # Aggregating this overall observations in "select" and divide by k
    # we get the formula R11 - ... 
    # R11 is the aggregation of the term "1" and divide by k
    Rint1 <- as.numeric((R11- t(((n-x_rank+1)/k)^gamma_1) %*% select/k)/gamma_1)

    # Rint2 is the estimator of the integral from 0 to 1 on R(1,t)/t
    # The rationale is symmetric as above, but take y dimension.
    # In addition, the integral of 1/t is log(t)    
    Rint2 <- as.numeric(- t(log((n-y_rank+1)/k)) %*% select/k)

    # Calculate R_g
    R_g <- R11 - gamma_1/(gamma_1+1)*((2*gamma_1+1)*Rint1-Rint2)

    # Calcluate gamma_tilde and variance
    gamma_tilde <- gamma_1-(gamma_1+1)*R_g*gamma_2
    var_gammat <- (1+gamma_tilde)^2*(1-m*R_g^2/(n+m))

    # Rint11 is an estimator of the integral from 0 to 1 on R(s,1)/s^{1-gamma}
    # The formula is the same as Rint1, with taking gamma=gamma_tilde
    # Since we have gamma_tilde, we can estimate it better
    Rint11 <- as.numeric((R11- t(((n-x_rank+1)/k)^gamma_tilde) %*% select/k)/gamma_tilde)

    # Rint3 is an estimator of the integral from 0 to 1 on R(1,t)/t*log(t)
    # The rationale is the same as Rint2, focusing on the y dimension.
    # In addition, the integral of 1/t* (log (t) is (log(t))^2/2
    Rint3 <- as.numeric(- t((log((n-y_rank+1)/k))^2) %*% select/k/2)


    # Calculate R_s, which is the term hat S_g/2 in the paper
    # It is used to calcluate the variance of hat sigma_g (with g=0)
    # Here we use gamma_tilde instead of gamma_1 since it is a better estimator
    R_s=(((3*gamma_tilde)-1)*Rint2+gamma_tilde*Rint3-((2*gamma_tilde)+1)^2*Rint11+2*(gamma_tilde+2)*R11)/2

    # Calculate sigma_tilde
    sigma_tilde <-sigma_1*(1-R_s*(sigma_2-1))

    return(list(gamma_origin = gamma_1, sigma_origin = sigma_1, gamma_semi = gamma_tilde, sigma_semi = sigma_tilde))
}

# Functions to calculate the variance of a quantile estimator
# For the traditional MLE


q_gamma_fun <- function (gamma, d){
        return(1/gamma*d^gamma*log(d)-(d^gamma-1)/(gamma^2))
}

var_Q_fun <- function (gamma){
    if(gamma>=0){
        return((gamma+1)^2)
    } else{
        return(1+4*gamma+5*gamma^2+2*gamma^3)
    }
}


# Improved quantile, with g=0
# data should be a two-column dataframe
# the first column, first n row, provides the main variable to be considered
# the second column, provides the covariate with more observations
# the first n cases match as pairs
# p is the tail probability level
# n is the sample size for the main variable
# k is a single value indicating the value k to be used


imp_quantile_mle <- function(data,p,n,k,k1){
    # Obtain m
    m <- dim(data)[1] - n
    
    # Obtain x
    x <- as.vector(data[1:n,1])
    x_rank <- rank(x)

    # Prepare y_tilde
    y_tilde <- tilde(as.vector(data[,2]),n)
    y_rank <- rank(y_tilde)
    y_origin <-data[,2]
    

    # Get the thresholds
    x_threshold=sort(x)[n-k]
    y_threshold=sort(y_tilde)[n-k]
    y_origin_threshold=sort(y_origin)[m+n-k1]

    # MLE
    ml <- gp.fit(x,x_threshold,method="zhang") 
    gamma_1 <- ml$approx.mean[2] 
    sigma_1 <- ml$approx.mean[1] 
    

        
    ml1 <- gp.fit(y_tilde,y_threshold,method="zhang")
    gamma_2 <- ml1$approx.mean[2]
    sigma_2 <- ml1$approx.mean[1]

    ml2 <- gp.fit(y_origin,y_origin_threshold,method="zhang") 
    gamma_origin <- ml2$approx.mean[2] 
    sigma_origin <- ml2$approx.mean[1] 
    var_gamma_origin <- (1+gamma_origin)^2

    # Calculate R11 and three integrals, see imp_gamma_mle for more explanation
    select <- (x>x_threshold & y_tilde>y_threshold)
    R11 <- sum(select)/k
    Rint1 <- as.numeric((R11- t(((n-x_rank+1)/k)^gamma_1) %*% select/k)/gamma_1)
    Rint2 <- as.numeric(- t(log((n-y_rank+1)/k)) %*% select/k)
    Rint3 <- as.numeric(- t((log((n-y_rank+1)/k))^2) %*% select/k/2)

    # Calculate R_g
    R_g <- R11 - gamma_1/(gamma_1+1)*((2*gamma_1+1)*Rint1-Rint2)

    # Calcluate gamma_tilde and variance
    gamma_tilde <- gamma_1-(gamma_1+1)*R_g*gamma_2
    var_gammat <- min((1+gamma_1)^2,(1+gamma_tilde)^2)*(1-m*R_g^2/(n+m))
    var_gamma <- (1+gamma_1)^2

    # Calculate R_s, which is the term hat S_g/2 in the paper
    Rint11 <- as.numeric((R11- t(((n-x_rank+1)/k)^gamma_tilde) %*% select/k)/gamma_tilde)
    R_s=(((3*gamma_tilde)-1)*Rint2+gamma_tilde*Rint3-((2*gamma_tilde)+1)^2*Rint11+2*(gamma_tilde+2)*R11)/2

    # Calculate sigma_tilde
    sigma_tilde <-sigma_1*(1-R_s*(sigma_2-1))

    # Calculate the quantiles
    # Original estimator
    Q <- x_threshold + sigma_1*((k/(n*p))^gamma_1-1)/gamma_1
    # Imprived estimator
    Q_tilde <- x_threshold + sigma_tilde*((k/(n*p))^gamma_tilde-1)/gamma_tilde
    # Original estimator for y
    Q_origin <- y_origin_threshold + sigma_origin*((k1/((n+m)*p))^gamma_origin-1)/gamma_origin

    # Calculate the variance of the quantile, for X and Y
    d=k/(n*p)
    d1=k1/((n+m)*p)
    var_Q=var_Q_fun(gamma_1)
    var_Q=var_Q/k*sigma_1^2*(q_gamma_fun(gamma_1,d))^2

    var_Q_origin = var_Q_fun(gamma_origin)
    var_Q_origin = var_Q_origin/k1 *sigma_origin^2*(q_gamma_fun(gamma_origin,d1))^2

    # Calculate the variance of the quantile for the semisupervised estimator
    # Recall that R_s is the term hat S_g/2 in the paper
    if(gamma_tilde>=0){
        var_Q_tilde = min((1+gamma_1)^2,(1+gamma_tilde)^2)*(1-m*R_g^2/(n+m))
    } else {
        # Calculate the terms Q1 and Q2 in the paper
        # Qint1 is Q1*(1+gamma), which simplifies also the later calcuation in Q
        Qint1 = -(gamma_tilde+2)*R11-gamma_tilde*Rint3+(-2*gamma_tilde+1)*Rint2+(2*gamma_tilde+1)*gamma_tilde*Rint1
        Qint2 = -2*R11+(2*gamma_tilde+1)^2/(gamma_tilde+1)*Rint1-gamma_tilde/(gamma_tilde+1)*Rint2
        # Calculate the term Q in the paper
        Qint = R_g*R_s+R_g*Qint1+R_s*Qint2

        # Calculate the term M1 and M2 in the paper
        Mint1 = (1+gamma_tilde)*(2*Rint2+Rint3-R11)
        Mint2 = - R_s*(3*Rint2+Rint3-2*R11)

        var_Q_tilde = (1+gamma_tilde)^2*(1+2*gamma_tilde)-m/(n+m)*((1+gamma_tilde)^2*R_g^2+gamma_tilde^2*R_s^2*2-2*gamma_tilde*(1+gamma_tilde)*Qint-2*gamma_tilde^2*Mint1+2*gamma_tilde^3*Mint2)
    }
    var_Q_tilde = var_Q_tilde/k * min(sigma_1^2,sigma_tilde^2)*min((q_gamma_fun(gamma_1,d))^2,(q_gamma_fun(gamma_tilde,d))^2)

    

    res_x_org = list(gamma=gamma_1, var_gamma=var_gamma, Q=Q, var_Q=var_Q)
    res_x_semi = list(gamma=gamma_tilde, var_gamma=var_gammat, Q=Q_tilde, var_Q=var_Q_tilde)
    res_y_org = list(gamma=gamma_origin, var_gamma=var_gamma_origin, Q=Q_origin, var_Q=var_Q_origin)
    return(list(res_x_org = res_x_org, res_x_semi = res_x_semi, res_y_org= res_y_org, R11=R11))
}
