Q <- length(tau)

jags_2mqrjm_n.weib.value <- function (...) {
  for(q in 1:Q){
    # constants
    c1[q] <- (1-2*tau[q])/(tau[q]*(1-tau[q]))
    c2[q] <- 2/(tau[q]*(1-tau[q]))
  }

  # likelihood
  for (i in 1:I){

    for(q in 1:Q){
    # longitudinal part
    for(j in offset[i]:(offset[i+1]-1)){
      # define object
      W[q, j] ~ dexp(1/sigma[q])
      prec[q,j] <- 1/(W[q, j]*sigma[q]*c2[q])
      # quantile distribution for quantile q
      y[q, j] ~ dnorm(mu[q, j], prec[q, j])
      mu[q, j] <- inprod(beta[q, 1:ncX], X[j, 1:ncX]) + inprod(b[q, i, 1:ncU], U[j, 1:ncU]) + c1[q]*W[q, j]
    }#end of j loop
    # random effects
    b[q, i, 1:ncU] ~ dmnorm(mu0[], prec.Sigma2[q, , ])
    # survival part
    shareQ[q, i] <- inprod(beta[q, 1:ncX], Xtime[i, 1:ncX]) + inprod(b[q, i, 1:ncU], Utime[i, 1:ncU])
    }
    etaBaseline[i] <- inprod(alpha[1: ncZ], Z[i, 1:ncZ])
  #fin boucle sur quantiles
    log_h1[i] <- log(shape) + (shape - 1) * log(Time[i]) + etaBaseline[i] + alpha.assoc * (shareQ[Q, i] - shareQ[Q-1, i])
    for (k in 1:K) {
      for (q in 1:Q){
        shareQ.s[q, i, k] <- inprod(beta[q, 1:ncX], Xs[K * (i - 1) + k, 1:ncX]) + inprod(b[q, i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])
      }
        SurvLong[i, k] <- wk[k] * shape * pow(st[i, k], shape - 1) * exp(alpha.assoc * (shareQ.s[Q, i, k] - shareQ.s[Q-1, i, k]) )
    }

    log_S1[i] <- (-exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ]))
    logL[i] <- event[i]*log_h1[i] + log_S1[i]
    mlogL[i] <- -logL[i] + C
    zeros[i] ~ dpois(mlogL[i])
  }#end of i loop
  # priors for parameters
  for (q in 1:Q){
  prec.Sigma2[q, 1:ncU, 1:ncU] ~ dwish(priorR.Sigma2[, ], priorK.Sigma2)
  covariance.b[q] <- inverse(prec.Sigma2[q, , ])

  beta[q, 1:ncX] ~ dmnorm(priorMean.beta[q, ], priorTau.beta[, ])
  sigma[q] ~ invgamma::dinvgamma(1/priorA.sigma, 1/priorB.sigma)
  }
  # priors for survival parameters
  alpha[1:ncZ] ~ dmnorm(priorMean.alpha[], priorTau.alpha[, ])
  shape ~ dgamma(priorA.shape, priorB.shape)
  alpha.assoc ~ dnorm(priorMean.alpha.assoc, priorTau.alpha.assoc)
}
