#######################################Inter test#####################################
IndaGE.inter <- function(data, pow = c(1:6, Inf), n.perm =1000, 
                              mainweights=weights1, interweights=weights2){
  ##### IndaGE for interaction test 
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), 
                             data = tdat, family = "binomial", verbose = F, control = list(sigma=1)) #, control = list(sigma=1)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  phi <- as.numeric(varc[nrow(varc),1])
  tau <- as.numeric(varc[1, 1])
  nullmodel$phi <- phi
  nullmodel$tau <- tau
  ################################################### 
  beta1 <- as.vector(nullmodel$coefficients$random$id) # random effects
  beta0 <- as.vector(nullmodel$coefficients$fixed)     #  fixed effects
  tmp0 <- as.vector(exp(G.weight%*%beta1))
  
  id <- which(is.infinite(tmp0))
  if(length(id)>0){
    tmp0 <- tmp0[-id] 
    G <- G[-id, ]
    G.weight <- G.weight[-id, ]
  }
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  tmp1 <- as.vector(exp(Z%*%beta0)) 
  tmp2 <- 1+as.vector(exp(Z%*%beta0))*mean(tmp0)
  
  EGY <- tmp1%*%t(colMeans(G*tmp0))/tmp2
  S.star <- as.vector(t(EG)%*%Y)-colSums(E*EGY)*weights2
  S.star <- S.star/phi
  Sigma.star <- GetVar_inter(Y, G, X, nullmodel, weights1)
  Sigma.star <- diag(weights2)%*%Sigma.star%*%diag(weights2)
  ############################################### 
  IndaGE <- GetInter_pval(Sigma=Sigma.star, U=S.star, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE 
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
   
  return(p.IndaGE)
}

##############################################

IndrareGE.inter <- function(data, n.perm =1000, 
                            mainweights=weights1, interweights=weights2){
  ##### IndrareGE for interaction test 
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), 
                             data = tdat, family = "binomial", verbose = F, control = list(sigma=1)) #, control = list(sigma=1)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  phi <- as.numeric(varc[nrow(varc),1])
  tau <- as.numeric(varc[1, 1])
  nullmodel$phi <- phi
  nullmodel$tau <- tau
  ################################################### 
  beta1 <- as.vector(nullmodel$coefficients$random$id) # random effects
  beta0 <- as.vector(nullmodel$coefficients$fixed)     #  fixed effects
  tmp0 <- as.vector(exp(G.weight%*%beta1))
  
  id <- which(is.infinite(tmp0))
  if(length(id)>0){
    tmp0 <- tmp0[-id] 
    G <- G[-id, ]
    G.weight <- G.weight[-id, ]
  }
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  tmp1 <- as.vector(exp(Z%*%beta0)) 
  tmp2 <- 1+as.vector(exp(Z%*%beta0))*mean(tmp0)
  
  EGY <- tmp1%*%t(colMeans(G*tmp0))/tmp2
  S.star <- as.vector(t(EG)%*%Y)-colSums(E*EGY)*weights2
  S.star <- S.star/phi
  Sigma.star <- GetVar_inter(Y, G, X, nullmodel, weights1)
  Sigma.star <- diag(weights2)%*%Sigma.star%*%diag(weights2)
  ############################################### 
  p.IndrareGE <- GetrareGE.inter_pval(Sigma.star, S.star)
  return(p.IndrareGE)
}

##############################################
aGE.inter <- function(data, pow = c(1:6, Inf), n.perm =1000, 
                         mainweights=weights1, interweights=weights2){
 
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
   
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
   
  MAFs <- colMeans(G, na.rm = T)/2
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  G.weight <- t(t(G)*weights1)
  
  ######################################
  aGE <- aGE::aGE(Y, G.weight, cov=X, model = "binomial", pow = pow, n.perm = n.perm, nonparaE = F)
  p.aGE <- aGE[1:9]
  return(p.aGE)
}

##############################################
rareGE.inter <- function(data, n.perm =1000, 
                mainweights=weights1, interweights=weights2){
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
 
  p.rareGE <-  rareGE::INT_RAN(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                               interweights = interweights, family = "binomial", binomialimpute = T)

  return(p.rareGE)
}

##############################################
WaGE.inter <- function(data, p.aGE, p.IndaGE, n.perm = n.perm, alpha=0.05){
  
  if(length(p.aGE)!=length(p.IndaGE)){stop("p.aGE and p.IndaGE should have the same length!")}
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
 
  re <- GetInd_pval(Y=Y, G=G, X=X, alpha, n.perm = n.perm) # independent test between E and G
  rho <- re$rho; p0 <- re$p0
  
  if(rho==1){
    p.waGE <- p.aGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
  }
  names(p.waGE) <- c(paste("waGEsm", 1:(length(p.aGE)-2), sep = ""), "waGEsm","waGEsm_fisher") 
  return(p.waGE)
}

################################################

WrareGE.inter <- function(data, p.rareGE, p.IndrareGE, n.perm = n.perm,  alpha=0.05){
  
  if(length(p.rareGE)!=length(p.IndrareGE)){stop("p.aGE and p.IndaGE should have the same length!")}
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  MAFs <- colMeans(G)/2
 
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  
  if(rho==1){
    p.wrareGE <- p.rareGE
 
  }else{
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  return(p.wrareGE)
}

################################################

Get.phi <- function(data, mainweights=weights1, interweights=weights2){
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  EG <- t(t(G)*weights2)*E
  G.weight <- t(t(G)*weights1)
  
  tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))
  if (k == 1) {
    fixedmod <- "trait ~ E"
  } else {
    fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
  if (m == 1) {
    randommod <- "~ 0 + G"} else {
      randommod <- paste("~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
    }
  
  nullmodel <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), 
                             data = tdat, family = "binomial", verbose = F)#, control = list(sigma=1)
  varc <- tryCatch(nlme::VarCorr(nullmodel), error=function(c) NA, warning = function(c) NA)
  tau <- as.numeric(varc[1, 1])
  phi <- as.numeric(varc[nrow(varc),1])
  
  return(phi)
}

#############################################Joint test#######################################

IndaGE.joint <- function(data, pow = c(1:6, Inf), n.perm =1000, 
                         mainweights=weights1, interweights=weights2){
  ##### IndaGE for joint test 
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  G.weight <- t(t(G)*weights1)
  
  weights <- c(weights1, weights2)
  G.joint <- cbind(G, G*E)
  G.joint.weight <- t(t(G.joint)*weights)
  fit0 <- glm(Y ~ X,family = "binomial")
  mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  res <- Y-mu  
  Sigma.star <- diag(weights)%*%GetVar_joint(G, Z, E, mu)%*%diag(weights)
  S.star <- as.vector(t(G.joint.weight)%*%Y)-c(g*sum(mu), g*sum(E*mu))*weights
  ######################################
  IndaGE <- GetJoint_pval(Sigma=Sigma.star, U=S.star, X=X, res=res, G=G, EG=G.joint.weight, n=n, pow=pow, n.perm=n.perm)
  p.IndaGE <- IndaGE$pvs
  names(p.IndaGE) <- c(paste("IndaGEsm", pow, sep = ""), "IndaGEsm","IndaGEsm_fisher") 
  
  return(p.IndaGE)
}

################################################################################
IndrareGE.joint <- function(data, rho=seq(0, 1, 0.1), n.perm =1000, 
                         mainweights=weights1, interweights=weights2){
 #### IndrareGE for joint test
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  G.weight <- t(t(G)*weights1)
  
  weights <- c(weights1, weights2)
  G.joint <- cbind(G, G*E)
  G.joint.weight <- t(t(G.joint)*weights)
  fit0 <- glm(Y ~ X,family = "binomial")
  mu <- as.vector(fit0$fitted.values)   ## the mean of Y
  res <- Y-mu  
  Sigma.star <- diag(weights)%*%GetVar_joint(G, Z, E, mu)%*%diag(weights)
  S.star <- as.vector(t(G.joint.weight)%*%Y)-c(g*sum(mu), g*sum(E*mu))*weights
  ######################################
  IndrareGE <- GetrareGE.joint_pval2(Uz=S.star, V=Sigma.star, n=n, rho=seq(0, 1, 0.1), n.perm=n.perm)
  p.IndrareGE <- IndrareGE$pval
  
  return(p.IndrareGE)
}

################################################################################a
aGE.joint <- function(data, pow = c(1:6, Inf), n.perm =1000, 
                            mainweights=weights1, interweights=weights2){
  #### IndrareGE for joint test
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  MAFs <- colMeans(G)/2
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  Z <- cbind(1, X)
  G.weight <- t(t(G)*weights1)
  
  p.aGE <- aGE::aGE.joint(Y, G.weight, cov=X, model="binomial", pow=pow, n.perm=n.perm, method="Simulation")
  return(p.aGE)
}

################################################################################
rareGE.joint <- function(data, rho=seq(0, 1, 0.1), n.perm =1000, 
                            mainweights=weights1, interweights=weights2){
  #### IndrareGE for joint test
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  
  MAFs <- colMeans(G, na.rm = T)/2
  Sigma <- var(G)
  g <- MAFs*2
  
  n <- dim(G)[1]
  m <- dim(G)[2]
  k <- ncol(X)
  
  if(is.function(mainweights)){
    weights1 <- mainweights(MAFs)
  }else{
    if(!is.vector(mainweights)) stop("Check the class of your variable mainweights: should be function or vector!")
    if(length(mainweights)!=m) stop("Number of variants inconsistent between genotypes and mainweights. Check your data...")
    weights1 <- mainweights
  }
  if(is.function(interweights)) {
    weights2 <- interweights(MAFs)
  }else{
    if(!is.vector(interweights)) stop("Check the class of your variable interweights: should be function or vector!")
    if(length(interweights)!=m) stop("Number of variants inconsistent between genotypes and interweights. Check your data...")
    weights2 <- interweights
  }  
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  
  ################################################### 
  rareGE <-  rareGE::JOINT(phenotype=Y, genotypes=G, covariates=X, mainweights = mainweights, 
                           interweights = interweights, family = "binomial", binomialimpute = T, 
                           rho = seq(0, 1, by = 0.1), B = n.perm)
  p.rareGE <- rareGE$pJOINT
  return(p.rareGE)
}

##############################################
WaGE.joint <- function(data, p.aGE, p.IndaGE, n.perm = n.perm, alpha=0.05){
  
  if(length(p.aGE)!=length(p.IndaGE)){stop("p.aGE and p.IndaGE should have the same length!")}
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
 
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm) # independent test between E and G
  rho <- re$rho; p0 <- re$p0
  
  if(rho==1){
    p.waGE <- p.aGE
  }else{
    t1 <- (1-p0)*(p.aGE-p.IndaGE)^2
    k1 <- t1/(t1+p0*p.aGE*p.IndaGE)
    TwaGE <- k1*tan((0.5-p.aGE)*pi)+(1-k1)*tan((0.5-p.IndaGE)*pi)
    p.waGE <- pcauchy(as.numeric(TwaGE), scale = 1, lower.tail = FALSE)
  }
  names(p.waGE) <- c(paste("waGEsm", 1:(length(p.aGE)-2), sep = ""), "waGEsm","waGEsm_fisher") 
  return(p.waGE)
}

################################################

WrareGE.joint <- function(data, p.rareGE, p.IndrareGE, n.perm = n.perm,  alpha=0.05){
  
  if(length(p.rareGE)!=length(p.IndrareGE)){stop("p.aGE and p.IndaGE should have the same length!")}
  
  Y <- data$Y
  X <- as.matrix(data$X)  ## the first column is environmental factor to be tested
  G <- as.matrix(data$G)
  MAFs <- colMeans(G)/2
  for(i in 1:ncol(G)){
    if(is.na(MAFs[i])) G[is.na(G[, i]), i] <- rbinom(length(G[is.na(G[, i]), i]), 2, mean(G[, i], na.rm = T)/2)
  }
  
  MAFs <- colMeans(G)/2
  id <- which(MAFs==0)
  if(length(id)>0) G <- G[, -id] 
  MAFs <- colMeans(G)/2
  
  E <- as.vector(scale(X[, 1], center = T, scale = T))
  X <- cbind(E, X[, -1])
  
  re <- GetInd_pval(Y, G, X, alpha, n.perm = n.perm)
  rho <- re$rho; p0 <- re$p0
  
  if(rho==1){
    p.wrareGE <- p.rareGE
    
  }else{
    t2 <- (1-p0)*(p.rareGE-p.IndrareGE)^2
    k2 <- t2/(t2+p0*p.rareGE*p.IndrareGE)
    TwrareGE <- k2*tan((0.5-p.rareGE)*pi)+(1-k2)*tan((0.5-p.IndrareGE)*pi)
    p.wrareGE <- pcauchy(as.numeric(TwrareGE), scale = 1, lower.tail = FALSE)
  }
  return(p.wrareGE)
}
