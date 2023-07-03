## testing tie-oriented modeling (methods) ##

# loading data
data(tie_reh)

# specifying linear predictor
tie_model <- ~ 1 + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(tie_mle)
expect_inherits(tie_mle,"remstimate")  
expect_length(tie_mle,18)
expect_identical(names(tie_mle),c("coefficients","loglik","gradient","hessian","vcov","se","diagnostics","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(tie_mle),8)
expect_identical(names(attributes(tie_mle)),c("names","class","formula","model","ordinal","method","approach","statistics"))
expect_identical(attr(tie_mle,"approach"),"Frequentist")
expect_silent(print(tie_mle))
expect_silent(summary(tie_mle))
expect_silent(residuals(tie_mle))
expect_silent(predict(tie_mle))
expect_silent(plot(tie_mle,tie_reh))
expect_silent(aic(tie_mle))
expect_silent(aicc(tie_mle))
expect_silent(bic(tie_mle))
expect_silent(waic(tie_mle))


# (2) method = "GDADAMAX" 
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10L))
tie_gdadamax <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10L)
expect_silent(tie_gdadamax)
expect_inherits(tie_gdadamax,"remstimate")  
expect_length(tie_gdadamax,18)
expect_identical(names(tie_gdadamax),c("coefficients","loglik","gradient","hessian","vcov","se","diagnostics","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(tie_gdadamax),11)
expect_identical(names(attributes(tie_gdadamax)),c("names","class","formula","model","ordinal","method","approach","statistics","epochs","epsilon","init"))
expect_identical(attr(tie_gdadamax,"approach"),"Frequentist")
expect_silent(print(tie_gdadamax))
expect_silent(summary(tie_gdadamax))
expect_silent(residuals(tie_gdadamax))
expect_silent(predict(tie_gdadamax))
expect_silent(plot(tie_gdadamax,tie_reh))
expect_silent(aic(tie_gdadamax))
expect_silent(aicc(tie_gdadamax))
expect_silent(bic(tie_gdadamax))
expect_silent(waic(tie_gdadamax))                   

# (3) method  = "BSIR" 

## (3.1) with prior = NULL
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL))
tie_bsir_no_prior <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL)
expect_silent(tie_bsir_no_prior)
expect_inherits(tie_bsir_no_prior,"remstimate")  
expect_length(tie_bsir_no_prior,13)
expect_identical(names(tie_bsir_no_prior),c("log_posterior","draws","log_proposal","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual","diagnostics"))
expect_length(attributes(tie_bsir_no_prior),10)
expect_identical(names(attributes(tie_bsir_no_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed"))
expect_identical(attr(tie_bsir_no_prior,"approach"),"Bayesian")
expect_silent(print(tie_bsir_no_prior))
expect_silent(summary(tie_bsir_no_prior))
expect_silent(residuals(tie_bsir_no_prior))
expect_silent(predict(tie_bsir_no_prior))
expect_silent(plot(tie_bsir_no_prior,tie_reh))
expect_error(aic(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(waic(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE) 

## (3.2) with a specified prior 
priormvt <- mvnfast::dmvt
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = priormvt,
                        mu = rep(0,dim(tie_reh_stats)[3]),
                        sigma = diag(dim(tie_reh_stats)[3]),
                        df = 1,
                        log = TRUE
                        ))
tie_bsir_with_prior <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = priormvt,
                        mu = rep(0,dim(tie_reh_stats)[3]),
                        sigma = diag(dim(tie_reh_stats)[3]),
                        df = 1,
                        log = TRUE)

expect_silent(tie_bsir_with_prior)
expect_inherits(tie_bsir_with_prior,"remstimate")  
expect_length(tie_bsir_with_prior,14)
expect_identical(names(tie_bsir_with_prior),c("log_posterior","draws","log_proposal","log_prior","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual","diagnostics"))
expect_length(attributes(tie_bsir_with_prior),11)
expect_identical(names(attributes(tie_bsir_with_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","prior","nsim","seed"))
expect_identical(attr(tie_bsir_with_prior,"approach"),"Bayesian")
expect_silent(print(tie_bsir_with_prior))
expect_silent(summary(tie_bsir_with_prior))
expect_silent(residuals(tie_bsir_with_prior))
expect_silent(predict(tie_bsir_with_prior))
expect_silent(plot(tie_bsir_with_prior,tie_reh))
expect_error(aic(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(waic(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE) 

# (4) method  = "HMC"    
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L)) 
tie_hmc <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L)
expect_silent(tie_hmc)
expect_inherits(tie_hmc,"remstimate")  
expect_length(tie_hmc,11)
expect_identical(names(tie_hmc),c("draws","log_posterior","coefficients","post.mean","vcov","sd","loglik","df.null","df.model","df.residual","diagnostics"))
expect_length(attributes(tie_hmc),14)
expect_identical(names(attributes(tie_hmc)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed","nchains","burnin","thin","init"))
expect_identical(attr(tie_hmc,"approach"),"Bayesian")
expect_silent(print(tie_hmc))
expect_silent(summary(tie_hmc))
expect_silent(residuals(tie_hmc))
expect_silent(predict(tie_hmc))
expect_silent(plot(tie_hmc,tie_reh))
expect_error(aic(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(waic(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE) 


# ordinal likelihood (tie-oriented modeling)
attr(tie_reh,"ordinal") <- TRUE
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))              
ordinal_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(print(ordinal_mle))
expect_silent(summary(ordinal_mle))  
expect_silent(residuals(ordinal_mle))
expect_silent(plot(ordinal_mle,tie_reh))


# testing estimation methods with active riskset 
tie_reh <- remify::remify(edgelist = tie_reh$edgelist, 
                            model = "tie", 
                            riskset="active")
# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
# (2) method = "GDADAMAX" 
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10L))
# (3) method  = "BSIR" 

## (3.1) with prior = NULL
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL))     

## (3.2) with a specified prior 
priormvt <- mvnfast::dmvt
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = priormvt,
                        mu = rep(0,dim(tie_reh_stats)[3]),
                        sigma = diag(dim(tie_reh_stats)[3]),
                        df = 1,
                        log = TRUE
                        ))
# (4) method  = "HMC"    
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L)) 