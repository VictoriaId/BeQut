load("C:/Users/vi/Desktop/Seafile/Seafile/2024_Victoria Idier/Progress/Progress.Rdata")
ech <- progress[which((progress$ethnie=="Non_Asian")&(progress$rand=="Placebo")),]
id <- unique(ech$ID)
id[100]
ech_progress <- ech[which(ech$ID<=200),]
length(unique(ech_progress$ID))

library(BeQut)
library(lqmm)

test_mqrjm <- mqrjm(formFixed = pressure ~ obstime + age + sex,
                    formRandom = ~ obstime,
                    formGroup = ~ ID,
                    formSurv = Surv(time, event) ~ age + sex,
                    survMod = "weibull",
                    param = "value",
                    timeVar = "obstime",
                    data = ech_progress,
                    tau = c(0.5, 0.75),
                    RE_ind = TRUE,
                    n.chains = 3,
                    n.iter = 5000,
                    n.burnin = 500,
                    n.thin = 1,
                    n.adapt = 10000,
                    precision = 10,
                    C = 1000,
                    save_jagsUI = TRUE,
                    save_va = FALSE,
                    parallel = FALSE)


formFixed = pressure ~ obstime + age + sex
formRandom = ~ obstime
formGroup = ~ ID
formSurv = Surv(time, event) ~ age + sex
survMod = "weibull"
param = "value"
timeVar = "obstime"
data = ech_progress
tau = c(0.5, 0.75)
RE_ind = TRUE
n.chains = 3
n.iter = 5000
n.burnin = 500
n.thin = 1
n.adapt = 10000
precision = 10
C = 1000
save_jagsUI = TRUE
save_va = FALSE
parallel = FALSE

jagsUI::traceplot(out_jags, parameters = "beta")
jagsUI::traceplot(out_jags, parameters = "alpha")
jagsUI::traceplot(out_jags, parameters = "alpha.assoc")
jagsUI::traceplot(out_jags, parameters = "sigma")
jagsUI::traceplot(out_jags, parameters = "b1")
jagsUI::traceplot(out_jags, parameters = "b2")
jagsUI::traceplot(out_jags, parameters = "covariance.b1")
jagsUI::traceplot(out_jags, parameters = "covariance.b2")
jagsUI::traceplot(out_jags, parameters = "shape")

plot(out_jags)
coda::cumuplot(out)
coda::gelman.plot(out_jags$sims.list$beta)
?gelman.plot
coda::gelman.plot(tmp_model)
