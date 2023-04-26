#tdat <- data.frame(trait=Y, X, G=G.weight, id=rep(1, n))

#if (k == 1) {
#  fixedmod <- "trait ~ E"
#} else {
#  fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
#if (m == 1) {
#  randommod <- "(~ 1|G)"
#  } else {
#    randommod <- paste0("(1|", names(tdat)[(k + 2):(1 + k + m)], ")", collapse = " + ")
#  }
# formula <- paste(fixedmod, "+", randommod)

#nullmodel <- lme4::glmer(formula = as.formula(formula), data=tdat, family = "binomial")
# TMB <- glmmTMB::glmmTMB(formula = as.formula(formula), data=tdat, family = "binomial", REML = T)

tdat <- data.frame(trait=Y, cbind(E, X[, -1]), G=G.weight, id=rep(1, n))

if (k == 1) {
  fixedmod <- "trait ~ E"
} else {
  fixedmod <- paste("trait ~", paste(names(tdat)[2:(k+1)], collapse = " + "))}
if (m == 1) {
  randommod <- " ~ 0 + G"} else {
    randommod <- paste(" ~ 0 +", paste(names(tdat)[(k + 2):(1 + k + m)], collapse = " + "))
  }

PQL <- MASS::glmmPQL(fixed = as.formula(fixedmod), random = list(id = nlme::pdIdent(as.formula(randommod))), data = tdat, family = "binomial", 
                     verbose = F)  
varc <- nlme::VarCorr(PQL)
phi <- as.numeric(varc[nrow(varc), 1])
tau <- as.numeric(varc[1, 1])
mu <- as.numeric(predict(PQL, data=tdat, type="response", level=1))
res <- Y-mu
#tmpdat <- data.frame(trait=rep(Y, each=m), X=envir, G=as.vector(t(G.weight)), id=rep(c(1:m), n))
#tmpdat <- tmpdat[order(tmpdat$id), ]
#tmpdat$X.V3 <- factor(tmpdat$X.V3)

# if (k == 1) {
#   formumod <- "trait ~ X"
# } else {
#   formumod <- paste("trait ~", paste(names(tmpdat)[2:(k+1)], collapse = " + "))
# }
#formumod <- paste(formumod, "+", "(", 1, "|", "G", ")" )
#TMB <- glmmTMB::glmmTMB(formula=as.formula(formumod), data=tmpdat, family = binomial(link = "logit"), REML =T )

#mu <- apply(matrix(as.numeric(fitted(TMB, data=tmpdat)), ncol = m, nrow = n, byrow = T), 1, mean)
#res <-  Y-mu
#tau <- summary(TMB)$varcor$cond$G[1] 
#phi <- sigma(TMB)^2
