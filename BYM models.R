##### Convergence plots
#mar =c(5, 3, 1, 2) + 0.1
par(mfrow=c(4,2))
plot(main ="", ylab = "Intercept", ridge03_bym_samples[,"beta0"], type = "l")
plot(main = "", xlab = "Intercept", density(ridge03_bym_samples[,"beta0"]))

plot(main = "", ylab = paste(names(data.clean[5])), ridge03_bym_samples[,"beta[1]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[5])), density(ridge03_bym_samples[,"beta[1]"]))

plot(main = "", ylab = paste(names(data.clean[6])), ridge03_bym_samples[,"beta[2]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[6])), density(ridge03_bym_samples[,"beta[2]"]))

plot(main = "", ylab = paste(names(data.clean[7])), ridge03_bym_samples[,"beta[3]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[7])), density(ridge03_bym_samples[,"beta[3]"]))

plot(main = "", ylab = paste(names(data.clean[8])),  ridge03_bym_samples[,"beta[4]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[8])), density(ridge03_bym_samples[,"beta[4]"]))

plot(main = "", ylab = paste(names(data.clean[9])), ridge03_bym_samples[,"beta[5]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[9])), density(ridge03_bym_samples[,"beta[5]"]))

plot(main = "", ylab = paste(names(data.clean[10])), ridge03_bym_samples[,"beta[6]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[10])), density(ridge03_bym_samples[,"beta[6]"]))

plot(main = "", ylab = paste(names(data.clean[11])), ridge03_bym_samples[,"beta[7]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[11])), density(ridge03_bym_samples[,"beta[7]"]))

plot(main = "", ylab = paste(names(data.clean[12])), ridge03_bym_samples[,"beta[8]"] , type = "l")
plot(main = "", xlab = paste(names(data.clean[12])), density(ridge03_bym_samples[,"beta[8]"]))

plot(main ="", ylab = expression(sigma^2), ridge03_bym_samples[,"sigma"], type ="l")
plot(main ="", xlab = expression(sigma^2), density(ridge03_bym_samples[,"sigma"]))

plot(main ="", ylab = expression(lambda), ridge03_bym_samples[,"lambda"], type = "l")
plot(main ="", xlab = expression(lambda), density(ridge03_bym_samples[,"lambda"]))


############## Model Fit
par(mfrow=c(3,2), mar=c(5, 3, 1, 2) + 0.1)
plot(x = fitted_ridge01_df$Fitted.Cases ,y = fitted_ridge01_df$Observed.Cases,
     main = "CAR Models", xlab = "Fitted values", ylab = "Observed Cases (adjacency-based)"); abline(a=0,b=1,lwd=2)

plot(x = fitted_ridge01_bym_df$Fitted.Cases ,y = fitted_ridge01_bym_df$Observed.Cases,
     main = "BYM Models", xlab = "Fitted values", ylab = "Observed Cases"); abline(a=0,b=1,lwd=2)

plot(x = fitted_ridge02_df$Fitted.Cases ,y = fitted_ridge02_df$Observed.Cases,
     main = "", xlab = "Fitted values", ylab = "Observed Cases (k = 3)"); abline(a=0,b=1,lwd=2)
plot(x = fitted_ridge02_bym_df$Fitted.Cases ,y = fitted_ridge02_bym_df$Observed.Cases,
     main = "", xlab = "Fitted values", ylab = "Observed Cases"); abline(a=0,b=1,lwd=2)

plot(x = fitted_ridge03_df$Fitted.Cases ,y = fitted_ridge03_df$Observed.Cases,
     main = "", xlab = "Fitted values", ylab = "Observed Cases (k = 5)"); abline(a=0,b=1,lwd=2)
plot(x = fitted_ridge03_bym_df$Fitted.Cases ,y = fitted_ridge03_bym_df$Observed.Cases,
     main = "", xlab = "Fitted values", ylab = "Observed Cases"); abline(a=0,b=1,lwd=2)

#########################################
#    Ridge Regression using BYM        # 
#########################################

# Defining the usual fixed priors model ------------------------------------------------------
ridge_code_bym <- nimbleCode({ 
  
  # Defining our priors
  beta0 ~ dnorm(0, sd = 200)
  
  for(j in 1:p){
    beta[j] ~ dnorm(0, sd=1*sqrt(1/lambda))  # Betas for the covariates
  }
  
  # Half-Cauchy priors for conditional std-dev & tuning param (lambda)
  #nu ~ T(dt(0, 1, 1),0,) 
  sigmaa ~ T(dt(0, 1, 1),0,)
  sigmab ~ T(dt(0, 1, 1),0,)
  lambda ~ T(dt(0, 1, 1),0,)
  
  # Defining transformed parameters
  tau2 <- 1/(sigmaa^2)
  
  # CAR prior
  s[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau2, zero_mean = 1)
  
  # Defining likelihood
  for(i in 1:N) { 
    log(r[i]) <- beta0 + inprod(X[i, 1:p], beta[1:p]) + s[i] + v[i]
    mu[i] <- E[i]*r[i]
    y[i] ~ dpois(mu[i])
    v[i] ~ dnorm(0, sigmab)
  }
})

# Constants and initial values --------------------------------------------------------------------------

set.seed(1234)

# (1) Define constants that go into CAR prior into a list
Consts_ridge03_bym = list(N = nrow(data.clean),             # Number of areas
                          p = ncol(data.clean[5:ncol(data.clean)]),  # Number ofcovariates
                          X = as.matrix(data.clean[, c(5:ncol(data.clean))]),# Matrix of covariates
                          L = length(adj.5),                   # Edges (length of adj vector)
                          E = data.clean$Expected,
                          adj = adj.5,                         # Adjacency matrix
                          num = num.nb.5,                      # Number of neighbors per district
                          weights = rep(1, length(adj.5)))     # Giving weight of 1 for neighbors

# (2) Define dependent variable (outcome)
outcome <- list(y = data.clean$Cases) # Observed cases

# (3) Define initial values for the MCMC chain
Inits_ridge <- list(beta0=0, 
                    beta=rep(0, ncol(data.clean[5:ncol(data.clean)])),
                    sigmaa=1, 
                    sigmab=1, 
                    lambda=1, 
                    s = rnorm(nrow(data.clean)),
                    v = rnorm(nrow(data.clean)))

# rep(0, nrow(data.clean)))

# Defining model -----------------------------------------------------------------------------
# (1) Define ridge nimble model
ridge03_bym_model = nimbleModel(code=ridge_code_bym, 
                                constants=Consts_ridge03_bym, 
                                data=outcome,
                                inits=Inits_ridge)

# (2) Compile ridge model for nimble
ridge03_bym_Cmodel = compileNimble(ridge03_bym_model)

# (3) Configure the MCMC chain
ridge03_bym_conf <- configureMCMC(ridge03_bym_model, 
                                  monitors = c('beta0','beta', 'sigmaa', 'sigmab', 'lambda', 's', 'v'))
# (4) Do something
ridge03_bym_conf$printSamplers()

# (5) Build MCMC object
ridge03_bym_MCMC = buildMCMC(ridge03_bym_conf, enableWAIC=T)

# (6) Compile model with MCMC object
ridge03_bym_CMCMC = compileNimble(ridge03_bym_MCMC, project=ridge03_bym_Cmodel)

# (7) Run the model
ridge03_bym_runMCMC = runMCMC(ridge03_bym_CMCMC, nburnin=10000, niter=300000, nchains=2, thin=25, WAIC=T)

# (8) Check convergence
ridge03_bym_samples = rbind(ridge03_bym_runMCMC$samples[[1]], ridge03_bym_runMCMC$samples[[2]])

# (9) WAIC
ridge03_bym_runMCMC$WAIC

# (10) Convergence
# Beta0
par(mfrow=c(1,1), mar =c(5, 3, 1, 2) + 0.1)
plot(main ="Intercept", ridge03_bym_samples[,"beta0"], type = "l")

# Betas (Covariates) - log scale
par(mfrow=c(2,4), mar =c(5, 3, 1, 2) + 0.1)
plot(main = paste(names(data.clean[5])),  ridge03_bym_samples[,"beta[1]"] , type = "l") 
plot(main = paste(names(data.clean[6])),  ridge03_bym_samples[,"beta[2]"] , type = "l")
plot(main = paste(names(data.clean[7])),  ridge03_bym_samples[,"beta[3]"] , type = "l")
plot(main = paste(names(data.clean[8])),  ridge03_bym_samples[,"beta[4]"] , type = "l")
plot(main = paste(names(data.clean[9])),  ridge03_bym_samples[,"beta[5]"] , type = "l")
plot(main = paste(names(data.clean[10])), ridge03_bym_samples[,"beta[6]"] , type = "l")
plot(main = paste(names(data.clean[11])), ridge03_bym_samples[,"beta[7]"] , type = "l")
plot(main = paste(names(data.clean[12])), ridge03_bym_samples[,"beta[8]"] , type = "l")

# sigma squared for S
par(mfrow=c(2,1), mar=c(2,4.5,2.1,1)+0.1)
plot(main ="Sigma2 S", ridge03_bym_samples[,"sigmaa"], type ="l")
plot(main ="Sigma2 V", ridge03_bym_samples[,"sigmab"], type ="l")

# lambda
plot(main ="Lambda", ridge03_bym_samples[,"lambda"], type = "l")

###### Posterior densities

# Betas (Covariates) - log scale
par(mfrow=c(3,4), mar=c(2,4.5,2.1,1)+0.1)
plot(main ="Intercept", density(ridge03_bym_samples[,"beta0"]))
plot(main = paste(names(data.clean[5])),density(ridge03_bym_samples[,"beta[1]"])) 
plot(main = paste(names(data.clean[6])),density(ridge03_bym_samples[,"beta[2]"]))
plot(main = paste(names(data.clean[7])),density(ridge03_bym_samples[,"beta[3]"]))
plot(main = paste(names(data.clean[8])),density(ridge03_bym_samples[,"beta[4]"]))
plot(main = paste(names(data.clean[9])),density(ridge03_bym_samples[,"beta[5]"]))
plot(main = paste(names(data.clean[10])),density(ridge03_bym_samples[,"beta[6]"]))
plot(main = paste(names(data.clean[11])),density(ridge03_bym_samples[,"beta[7]"]))
plot(main = paste(names(data.clean[12])),density(ridge03_bym_samples[,"beta[8]"]))

#sigma squared for S
par(mfrow=c(1,2), mar=c(2,4.5,2.1,1)+0.1)
plot(main ="Sigma Squared S", density(ridge03_bym_samples[,"sigmaa"]))
plot(main ="Sigma Squared V", density(ridge03_bym_samples[,"sigmab"]))

#Lambda
plot(main ="Lambda", density(ridge03_bym_samples[,"lambda"]))

#Spatial effects for each region 
par(mfrow=c(4,5), mar=c(2,3,2.1,1)+0.1)
plot(main = paste(data.clean$District[1]), density(ridge03_bym_samples[,"s[1]"])) 
plot(main = paste(data.clean$District[2]), density(ridge03_bym_samples[,"s[2]"]))
plot(main = paste(data.clean$District[3]), density(ridge03_bym_samples[,"s[3]"]))
plot(main = paste(data.clean$District[4]), density(ridge03_bym_samples[,"s[4]"]))
plot(main = paste(data.clean$District[5]), density(ridge03_bym_samples[,"s[5]"]))
plot(main = paste(data.clean$District[6]), density(ridge03_bym_samples[,"s[6]"]))
plot(main = paste(data.clean$District[7]), density(ridge03_bym_samples[,"s[7]"]))
plot(main = paste(data.clean$District[8]), density(ridge03_bym_samples[,"s[8]"]))
plot(main = paste(data.clean$District[9]), density(ridge03_bym_samples[,"s[9]"]))
plot(main = paste(data.clean$District[10]), density(ridge03_bym_samples[,"s[10]"]))
plot(main = paste(data.clean$District[11]), density(ridge03_bym_samples[,"s[11]"]))
plot(main = paste(data.clean$District[12]), density(ridge03_bym_samples[,"s[12]"]))
plot(main = paste(data.clean$District[13]), density(ridge03_bym_samples[,"s[13]"]))
plot(main = paste(data.clean$District[14]), density(ridge03_bym_samples[,"s[14]"]))
plot(main = paste(data.clean$District[15]), density(ridge03_bym_samples[,"s[15]"]))
plot(main = paste(data.clean$District[16]), density(ridge03_bym_samples[,"s[16]"]))
plot(main = paste(data.clean$District[17]), density(ridge03_bym_samples[,"s[17]"]))
plot(main = paste(data.clean$District[18]), density(ridge03_bym_samples[,"s[18]"]))

# (11) Posterior summaries: Means and CrIs -------------------------------------------------------

betas <- paste0("beta[", 1:8, "]")
s <- paste0("s[", 1:18, "]")
v <- paste0("v[", 1:18, "]")

################## Beta0

mean_beta0_ridge03_bym <- mean(ridge03_bym_samples[ ,"beta0"])

#Calculating the Beta 95% CrI 
CrI_beta0_ridge03_bym <- t(quantile(ridge03_bym_samples[ , "beta0"], probs=c(0.025,0.975)))

#Combining the means and CIs into a dataframe
mean_CrI_beta0_ridge03_bym <- cbind(mean_beta0_ridge03_bym, CrI_beta0_ridge03_bym) %>% as.data.frame()

#Renaming the columns
colnames(mean_CrI_beta0_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Converting means and CrIs from log to regular scale
round(exp(mean_CrI_beta0_ridge03_bym),2)

################## Covariate Estimates: Beta

#Calculating the Beta means 
means_betas_ridge03_bym <- apply(ridge03_bym_samples[,1:8],2, mean)

#Calculating the Beta 95% CrI 
CrI_betas_ridge03_bym <- t(apply(ridge03_bym_samples[ ,1:8],2, function(x) quantile(x, probs=c(0.025,0.975))))

#Combining the means and CIs into a dataframe
mean_CrI_betas_ridge03_bym <- t(matrix(c(means_betas_ridge03_bym, CrI_betas_ridge03_bym), ncol=8,  byrow=T)) %>% as.data.frame()

#Renaming the columns
colnames(mean_CrI_betas_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Converting means and CrIs from log to regular scale
round(exp(mean_CrI_betas_ridge03_bym), 2)


################### Spatial effect per region: S
library(dplyr)
#Calculating the S means
means_s_ridge03_bym <- apply(ridge03_bym_samples[ ,s[1:18]],2, mean)

#Calculating the Beta 95% CrI 
CrI_s_ridge03_bym <- t(apply(ridge03_bym_samples[ , s[1:18]],2, 
                             function(x) quantile(x, probs=c(0.025,0.975))))

#Combining the means and CIs into a dataframe
mean_CrI_s_ridge03_bym <- cbind(means_s_ridge03_bym, CrI_s_ridge03_bym) %>% as.data.frame()

#Renaming the columns
colnames(mean_CrI_s_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Converting means and CrIs from log to regular scale
round(exp(mean_CrI_s_ridge03_bym),2)

################### Additional error variable: v

#Calculating the v means
means_v_ridge03_bym <- apply(ridge03_bym_samples[ ,v[1:18]],2, mean)

#Calculating the Beta 95% CrI 
CrI_v_ridge03_bym <- t(apply(ridge03_bym_samples[ , v[1:18]],2, 
                             function(x) quantile(x, probs=c(0.025,0.975))))

#Combining the means and CIs into a dataframe
mean_CrI_v_ridge03_bym <- cbind(means_v_ridge03_bym, CrI_v_ridge03_bym) %>% as.data.frame()

#Renaming the columns
colnames(mean_CrI_v_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Converting means and CrIs from log to regular scale
exp(mean_CrI_v_ridge03_bym)

################### Sigma a (1/tau)

#Calculating the sigma means
means_sigmaa_ridge03_bym <- mean(ridge03_bym_samples[ , "sigmaa"])

#Calculating the Beta 95% CrI 
CrI_sigmaa_ridge03_bym <- t(quantile(ridge03_bym_samples[ , "sigmaa"], probs=c(0.025,0.975)))

#Combining the means and CIs into a dataframe
mean_CrI_sigmaa_ridge03_bym <- cbind(means_sigmaa_ridge03_bym, CrI_sigmaa_ridge03_bym) %>%
  as.data.frame()

#Renaming the columns
colnames(mean_CrI_sigmaa_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Consigmaerting means and CrIs from log to regular scale
round(mean_CrI_sigmaa_ridge03_bym,2)

#Tau
mean_CrI_tau_ridge03_bym <- 1/(mean_CrI_sigmaa_ridge03_bym)
round(mean_CrI_tau_ridge03_bym,3)

################### Sigma b (independent error)

#Calculating the sigma means
means_sigmab_ridge03_bym <- mean(ridge03_bym_samples[ , "sigmab"])

#Calculating the Beta 95% CrI 
CrI_sigmab_ridge03_bym <- t(quantile(ridge03_bym_samples[ , "sigmab"], probs=c(0.025,0.975)))

#Combining the means and CIs into a dataframe
mean_CrI_sigmab_ridge03_bym <- cbind(means_sigmab_ridge03_bym, CrI_sigmab_ridge03_bym) %>%
  as.data.frame()

#Renaming the columns
colnames(mean_CrI_sigmab_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Consigmaerting means and CrIs from log to regular scale
round(mean_CrI_sigmab_ridge03_bym,2)

################## Lambda

mean_lambda_ridge03_bym <- mean(ridge03_bym_samples[ ,"lambda"])

#Calculating the Beta 95% CrI 
CrI_lambda_ridge03_bym <- t(quantile(ridge03_bym_samples[ , "lambda"], probs=c(0.025,0.975)))

#Combining the means and CIs into a dataframe
mean_CrI_lambda_ridge03_bym <- cbind(mean_lambda_ridge03_bym, CrI_lambda_ridge03_bym) %>% as.data.frame()

#Renaming the columns
colnames(mean_CrI_lambda_ridge03_bym) <- c("Mean", "95% CrI LL", "95% CrI UL")

# Converting means and CrIs from log to regular scale
exp(mean_CrI_lambda_ridge03_bym)

###### Model Fit
#Creating matrix to store fitted values
fitted_ridge03_bym = matrix(NA, nrow=nrow(ridge03_bym_samples), ncol= nrow(data.clean))

#Loop to generate fitted values: log(ri)
for(i in 1:nrow(data.clean)){
  fitted_ridge03_bym[,i] = ridge03_bym_samples[,"beta0"] + 
    as.matrix(data.clean[i,5:ncol(data.clean)]) %*%t(ridge03_bym_samples[,1:8]) + 
    ridge03_bym_samples[, s[i]]+ 
    ridge03_bym_samples[, v[i]] }

#Loop to generate fitted values: log(ri)
for(i in 1:nrow(data.clean)){
  fitted_ridge03_bym[,i] = ridge03_bym_samples[,"beta0"] + as.matrix(data.clean[i,5:ncol(data.clean)]) %*%t(ridge03_bym_samples[,1:8])+ ridge03_bym_samples[, s[i]]+ ridge03_bym_samples[, v[i]] }


# Mean fitted SIR
fitted_SIR_ridge03_bym <- exp(apply(fitted_ridge03_bym, 2, mean))

# Mean fitted Cases
fitted_ridge03_bym_cases <- fitted_SIR_ridge03_bym*data.clean$Expected

# 95% CrI SIR
fitted_ridge03_CrI_SIR <- t(exp(apply(fitted_ridge03_bym, 2, function(x) quantile(x, probs=c(0.025, 0.975)))))

# 95% CrI Cases
fitted_ridge03_bym_CrI <-  data.clean$Expected *fitted_ridge03_CrI_SIR

# Creating Dataframe to store fitted values and quantiles in regular svale
fitted_ridge03_bym_df <- data.frame("Observed SIR"   =  round(data.clean$SIR, 2), 
                                    "Fitted SIR"     =  round(fitted_SIR_ridge03_bym,2), 
                                    "Observed Cases" =  round(data.clean$Cases), 
                                    "Fitted Cases" =  round(fitted_ridge03_bym_cases), 
                                    "CrI LL"    =  round(fitted_ridge03_bym_CrI[,1]),
                                    "CrI UL"    =  round(fitted_ridge03_bym_CrI[,2]))


# Top 6 rows
head(fitted_ridge03_bym_df)









