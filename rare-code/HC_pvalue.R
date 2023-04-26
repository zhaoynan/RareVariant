
## To calculate the variance of r(t)
Get_bar.rk <- function(Sigma_Z, k){
  ## k is the order to sum ,k=1,2, 3,.......
  # Sigma_Z is the covriance of marginal score statistic Z
  
  Sigma_Z <- as.matrix(Sigma_Z)
  m <- dim(Sigma_Z)[1]
  
  Sigma_k <- Sigma_Z^k
  sum.r <- sum(Sigma_k)-sum(diag(Sigma_k))
  bar.rk <- sum.r/(m*(m-1))
  
  return(bar.rk)
}

## To obtain the Hermite matrix 
Get_Hermite <- function(x, k){
  ## k is the order to sum ,k=1,2,3,.......
  ## x  is a scalar or vector
  
  H <- matrix(NA, ncol = k+1, nrow = length(x))
  H[, 1] <- 1
  H[, 2] <- x
  
  for(i in 1:length(x)){
    for(j in (3:(k+1))){
      H[i,j] <- x[i]*H[i, (j-1)]-(j-2)*H[i, (j-2)]
    }
  }
  return(H)
}

# To estimate the vaiance of r(t)
Get_Var.rt <- function(Sigma_Z, t, k){
  ## k is the order to sum ,n=1,2, 3,.......
  ## t is the possible p value
  ## H is the Hermite matrix
  
  t[which(t==0)] <- 1e-16
  t[which(t==1)] <- 1-1e-16
  
  Sigma_Z <- as.matrix(Sigma_Z)
  m <- dim(Sigma_Z)[1]
  
  if(cor(as.vector(Sigma_Z),as.vector(diag(m)))>=0.95){
    var.rt <- m*t*(1-t)
  }else{
    if(any(t>1) || any(t<0)){
      print("t is between 0 and 1!!")
    }else{
      Z.t <- qnorm(t/2, mean=0, sd=1, lower.tail = FALSE)
      sum1 <- m*t*(1-t)
      sum2 <- 4*m*(m-1)*dnorm(Z.t, mean=0, sd=1)^2
      bar.rk <- fact <- numeric(k)
      for(i in 1:k){
        bar.rk[i] <- Get_bar.rk(Sigma_Z, 2*i)
        fact[i] <- factorial(2*i)
      }
      if(length(Z.t)==1){
        H.t <- matrix(Get_Hermite(Z.t, (2*k-1))[, seq(2, 2*k, 2)], nrow = 1)
      }else{
        H.t <- Get_Hermite(Z.t, (2*k-1))[, seq(2, 2*k, 2)]
      }
      
      sum3 <- colSums(t(H.t^2)*bar.rk/fact) 
      var.rt <- sum1+sum2*sum3
    }
  } 
  return(var.rt)
  
}

## To calculate the statistic of WRHC
stat_WRHC <- function(Sigma_Z, p_score, k){
  ## var.r is the variance of r
  ## weights is the weight vector
  
  p_score <- as.vector(p_score)
  if(any(p_score>1) || any(p_score<0)){
    print("p_score should be between 0 and 1!!")
  }else{
    m <- length(p_score)
    var.r <- Get_Var.rt(Sigma_Z, p_score, k)
    #r <- rank(p_score) 
    r <- numeric(m)
    for(i in 1:m){
      r[i] <- sum(p_score <= p_score[i])
    }
    
    
    HC.t <- (r-m*p_score)/sqrt(var.r)
    HC <- max(HC.t)
    
    return(HC)
  }
}


## To calculate tk 

Get_tk <- function(data, Sigma_Z, h, k){
  ## h is the calculated statistic
  ## weights is the weight vector
  G <- data$G
  
  m <- dim(Sigma_Z)[1]
  t <- numeric(m)
  
  #w.t <- weights[order(p_score)]
  
  for(i in 1:m){
    c.t <- function(t){
      var.rt <- Get_Var.rt(Sigma_Z, t, k)
      c.wrhc <- h*sqrt(var.rt)+t*m
      return(c.wrhc-i)
    }
    t[i] <- uniroot(c.t, c(1e-50, (1-1e-50)), tol =1e-50)$root
    if(t[i]==1) t[i] <- 1-1e-16
  }
  t <- sort(t)
  return(t) 
}

## beta-binomial distribution:betabin(s, alpha, beta)

Get_beta.bin <- function(x, s, alpha, beta){
  ## x is a non-negative integer
  ## s is a positive integer not less than x
  if(any(x>s)){
    print("s is a positive integer not less than x!")
  }else{
    comm <- factorial(s)/factorial(x)/factorial(s-x)
    Num <- beta(x+alpha, s-x+beta)
    Den <- beta(alpha, beta)
    if(Den==0){
      ret <- 0
    }else{
      ret <- comm*Num/Den
    }
    return(ret)
  }
}


## Calculate the probability P(|Z_j|, |Z_l| >= |Z_tk|)

Get_P.jl <- function(Sigma_Z, t, k){
  m <- dim(Sigma_Z)[1]
  Z.t <- qnorm(t/2, mean=0, sd=1, lower.tail = FALSE)
  
  if(length(Z.t)==1){
    H.t <- matrix(Get_Hermite(Z.t, (2*k-1))[, seq(2, 2*k, 2)], nrow = 1)
  }else{
    H.t <- Get_Hermite(Z.t, (2*k-1))[, seq(2, 2*k, 2)]
  }
  
  fact <- matrix(NA, nrow=m, ncol=k)
  ret <- list()
  
  for(i in 1:(m-1)){
    ret[[i]] <- matrix(NA, ncol=m, nrow=m-i)
    for(j in (i+1):m){
      fact[, 1:k] <- t((Sigma_Z[i, j]^(2*(1:k))/factorial(2*(1:k)))*t(H.t[, 1:k]^2))
      ret[[i]][(j-i), ] <- rowSums(fact)*dnorm(Z.t)^2+t^2/2
    }
  }
  return(ret)
}

## Estimate the parameter of beta-binomial(s, alpha, beta)

Get_par <- function(Sigma_Z, t, s, k){
  ## t is a vector with m
  # s is the value of t(k+1), s=1, 2, ....
  
  m <- dim(Sigma_Z)[1]
  
  var.rt <- Get_Var.rt(Sigma_Z, t, k)
  Pjl.tk <- Get_P.jl(Sigma_Z, t, k)
  E <- Var <- alpha <- beta <- list()
  
  for(i in 1:length(s)){
    
    E[[i]] <- Var[[i]] <- numeric(m)
    # the expectation and variance of r_t(k)|r_t(k+1)=s, k=1, 2,....., m-1
    E[[i]][-m] <- s[i]*t[-m]/t[-1]
    if(s[i]<=1){
      Var[[i]][-m] <- E[[i]][-m]-E[[i]][-m]^2
    }else{
      sum <- numeric(m-1)
      for(l in 1:(m-2)){
        sum <- sum+colSums(Pjl.tk[[l]][, -m]/Pjl.tk[[l]][, -1])
      }
      sum <- sum+Pjl.tk[[m-1]][, -m]/Pjl.tk[[m-1]][, -1]
      Var[[i]][-m] <- 2*(s[i]*(s[i]-1)/m/(m-1))*sum+E[[i]][-m]-E[[i]][-m]^2
    }
    
    E[[i]][m] <- m*t[m]
    Var[[i]][m] <- var.rt[m]
    
    alpha[[i]] <- beta[[i]] <- numeric(m)
    if(s[i]==0){
      alpha[[i]][-m] <- beta[[i]][-m] <- 0.5
    }
    if(s[i]==1){
      alpha[[i]][-m] <- 0.5
      beta[[i]][-m] <- (1-E[[i]][-m])*alpha[[i]][-m]/E[[i]][-m]
    }
    if(s[i]>1){
      h <- Var[[i]][-m]/(s[i]*(t[-m]/t[-1])*(1-t[-m]/t[-1]))
      #h <- Var[[i]][-m]*s[i]/E[[i]][-m]/(s[i]-E[[i]][-m])
      alpha[[i]][-m] <- (s[i]-h)*E[[i]][-m]/s[i]/(h-1)
      beta[[i]][-m] <- (s[i]/E[[i]][-m]-1)*alpha[[i]][-m]
    }
    #h <- Var[[i]][m]*m/E[[i]][m]/(m-E[[i]][m])
    h <- Var[[i]][m]/(m*t[m]*(1-t[m]))
    alpha[[i]][m] <- (m-h)*t[m]/(h-1)
    beta[[i]][m] <- (1/t[m]-1)*alpha[[i]][m]
    
  }
  return(list(alpha=alpha, beta=beta, E=E, Var=Var))  
}

Get_pvalue <- function(Sigma_Z, t, k){
  
  m <- dim(Sigma_Z)[1]
  q.ka <- matrix(NA, ncol=m, nrow=m)
  
  if(cor(as.vector(Sigma_Z),as.vector(diag(m)))>=0.95){
    q.ma <- dbinom(0:(m-1), m, t[m])
    q.ka[m, ] <- q.ma
    
    for(i in (m-1):1){
      pi <- matrix(0, ncol=i+1, nrow=i+1)
      for(s in 0:i){
        pi[s+1, 1:(s+1)] <- dbinom(0:s, s, t[i]/t[i+1])
        
      }
      q.ka[i, 1:i] <- colSums(as.matrix(pi[, -(i+1)]*q.ka[i+1, 1:(i+1)]), na.rm = TRUE)/sum(q.ka[i+1, ], na.rm=TRUE)
    }
  }else{
    par.betabin <- Get_par(Sigma_Z, t, m, k)
    alpha <- par.betabin$alpha[[1]]
    beta <- par.betabin$beta[[1]]
    
    #if(alpha[m]<0 || beta[m]<0){
    #  q.ma <- dbinom(0:(m-1), m, t[m])
    #  q.ka[m, ] <- q.ma
    #}else{
    q.ma <- Get_beta.bin(0:(m-1), m, alpha[m], beta[m])
    q.ka[m, ] <- q.ma
    #}
    
    for(i in (m-1):1){
      pi <- matrix(0, ncol=i+1, nrow=i+1)
      par.betabin <- Get_par(Sigma_Z, t, 0:i, k)
      for(s in 0:i){
        alpha <- par.betabin$alpha[[s+1]]
        beta <- par.betabin$beta[[s+1]]
        pi[s+1, 1:(s+1)] <- Get_beta.bin(0:s, s, alpha[i], beta[i])
      }
      q.ka[i, 1:i] <- colSums(as.matrix(pi[, -(i+1)]*q.ka[i+1, 1:(i+1)]), na.rm = TRUE)/sum(q.ka[i+1, ], na.rm=TRUE) 
    } 
  }
  sum.all <- as.matrix(rowSums(q.ka, na.rm=TRUE))
  p.value <- 1-apply(sum.all, 2, prod)
  
  return(list(q.ka=q.ka[!upper.tri(q.ka)], p.value=p.value))
}


