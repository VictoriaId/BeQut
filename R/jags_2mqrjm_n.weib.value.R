Q <- length(tau)

jags_2mqrjm_n.weib.value <- function (...) {
  for(q in 1:Q){
    # constants
    c1[q] <- (1-2*tau[q])/(tau[q]*(1-tau[q]))
    c2[q] <- 2/(tau[q]*(1-tau[q]))
  }

  # likelihood
  for (i in 1:I){

    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      # define object
      W[1, j] ~ dexp(1/sigma[1])
      W[2, j] ~ dexp(1/sigma[2])
      prec1[j] <- 1/(W[1, j]*sigma[1]*c2[1])
      prec2[j] <- 1/(W[2, j]*sigma[2]*c2[2])
      # quantile distribution for quantile q
      y[j, 1] ~ dnorm(mu1[j], prec1[j])
      mu1[j] <- inprod(beta[1, 1:ncX], X[j, 1:ncX]) + inprod(b1[i, 1:ncU], U[j, 1:ncU]) + c1[1]*W[1, j]
      y[j, 2] ~ dnorm(mu2[j], prec2[j])
      mu2[j] <- inprod(beta[2, 1:ncX], X[j, 1:ncX]) + inprod(b2[i, 1:ncU], U[j, 1:ncU]) + c1[2]*W[2, j]
    }#end of j loop
    # random effects
    b1[i, 1:ncU] ~ dmnorm(mu0[], prec1.Sigma2[ , ])
    b2[i, 1:ncU] ~ dmnorm(mu0[], prec2.Sigma2[ , ])
    # survival part
    shareQ1[i] <- inprod(beta[1, 1:ncX], Xtime[i, 1:ncX]) + inprod(b1[i, 1:ncU], Utime[i, 1:ncU])
    shareQ2[i] <- inprod(beta[2, 1:ncX], Xtime[i, 1:ncX]) + inprod(b2[i, 1:ncU], Utime[i, 1:ncU])
    etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
    log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alpha.assoc * (shareQ2[i] - shareQ1[i])
    for (k in 1:K) {
        shareQ1.s[i, k] <- inprod(beta[1, 1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b1[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
        shareQ2.s[i, k] <- inprod(beta[2, 1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b2[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
        SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alpha.assoc * (shareQ2.s[i, k] - shareQ1.s[i, k]) )
    }

    log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for parameters
  prec1.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b1 <- inverse(prec1.Sigma2[ , ])
  prec2.Sigma2[1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b2 <- inverse(prec2.Sigma2[ , ])

  for (q in 1:Q){
  beta[q, 1:ncX] ~ dmnorm(priorMean.beta[q, ], priorTau.beta[, ])
  invsigma[q] ~ dgamma(priorA.invsigma, priorB.invsigma)
  sigma[q] <- 1/invsigma[q]
  }
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  alpha.assoc ~ dnorm(priorMean.alpha.assoc, priorTau.alpha.assoc)
}
