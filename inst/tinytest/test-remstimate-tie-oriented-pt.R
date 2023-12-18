## testing tie-oriented modeling (methods) ##

# loading data
data(tie_data)

set_seed <- 23929

# processing data with simultaneous events (two or more events observed at the same time point)
tie_data$edgelist$time <- floor(tie_data$edgelist$time) 
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")

# specifying linear predictor
tie_model <- ~ 1  + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 

# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model, method = "pt")

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
# testing method = NULL (default is "MLE")
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = NULL))
                        
tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(tie_mle)
expect_inherits(tie_mle,"remstimate")  
expect_length(tie_mle,17)
expect_identical(names(tie_mle),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(tie_mle),10)
expect_identical(names(attributes(tie_mle)),c("names","class","formula","model","ordinal","method","approach","statistics","where_is_baseline","ncores"))
expect_identical(attr(tie_mle,"approach"),"Frequentist")
expect_silent(print(tie_mle))
expect_silent(summary(tie_mle))
expect_silent(diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_silent(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics, which = 2, effects = "inertia")) # plotting a specific effect
expect_silent(plot(x = tie_mle,reh = tie_reh, stats = tie_reh_stats)) # without supplying diagnostics but supplying the array of stats
expect_silent(aic(tie_mle))
expect_silent(aicc(tie_mle))
expect_silent(bic(tie_mle))
expect_silent(waic(tie_mle))

# error when 'diagnostics' has at least one statistic's name wrong
colnames(tie_reh_diagnostics$residuals$smoothing_weights)[1] <- "INDEGREEsender"
expect_error(plot(x = tie_mle,reh = tie_reh, diagnostics = tie_reh_diagnostics),
"one or more effects not found inside the object 'diagnostics'.",
fixed=TRUE)

# test diagnostics with only baseline
tie_stats_baseline <- remstats::remstats(reh = tie_reh, tie_effects = ~1 , method = "pt")
tie_mle_only_baseline <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_stats_baseline,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = tie_mle_only_baseline, reh = tie_reh, stats = tie_stats_baseline))


# tests on "WAIC" for "MLE"

# WAIC = TRUE
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE))
# WAIC = NULL
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = NULL))
# WAIC not a logical scalar
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = 20))

# nsimWAIC is not a numeric scalar                       
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = "text"))
# nsimWAIC is supplied                        
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
tie_mle_with_waic <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100)
expect_silent(print(tie_mle_with_waic))
expect_silent(summary(tie_mle_with_waic))    
expect_silent(waic(tie_mle_with_waic))                  


# ordinal likelihood: WAIC for ordinal likelihod + testing omit_dyad routine
omit_dyad <- list()
omit_dyad[[1]] <- list(time = c(120,158), dyad = data.frame(actor1=c(NA,"4"),actor2=c("4",NA),type=c(NA,NA))) # excluding actor 4 from the risk set between (observed) time points 120 and 158
tie_reh_ordinal <- remify::remify(edgelist = tie_data$edgelist, model = "tie", ordinal = TRUE, omit_dyad = omit_dyad)
tie_reh_stats_ordinal <- remstats::remstats(reh = tie_reh_ordinal, tie_effects = tie_model, method = "pt")
expect_silent(remstimate::remstimate(reh = tie_reh_ordinal,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
# interval likelihood: WAIC for ordinal likelihod + testing omit_dyad routine in estimation (remstimate) and diagnostics method
tie_reh_omit_test <- remify::remify(edgelist = tie_data$edgelist, model = "tie", omit_dyad = omit_dyad)
expect_silent(remstimate::remstimate(reh = tie_reh_omit_test,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
tie_mle_omit_test <- remstimate::remstimate(reh = tie_reh_omit_test,
                        stats = tie_reh_stats_ordinal,
                        ncores = 1L)                     
expect_silent(diagnostics(object = tie_mle_omit_test, reh = tie_reh_omit_test, stats = tie_reh_stats_ordinal))   
                     
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
expect_length(tie_gdadamax,17)
expect_identical(names(tie_gdadamax),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(tie_gdadamax),13)
expect_identical(names(attributes(tie_gdadamax)),c("names","class","formula","model","ordinal","method","approach","statistics","epochs","epsilon","init","where_is_baseline","ncores"))
expect_identical(attr(tie_gdadamax,"approach"),"Frequentist")
expect_silent(print(tie_gdadamax))
expect_silent(summary(tie_gdadamax))
expect_silent(diagnostics(object = tie_gdadamax, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_gdadamax, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_gdadamax,reh = tie_reh, diagnostics = tie_reh_diagnostics))
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
                        prior = NULL,
                        seed = set_seed))
tie_bsir_no_prior <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL,
                        seed = set_seed)
expect_silent(tie_bsir_no_prior)
expect_inherits(tie_bsir_no_prior,"remstimate")  
expect_length(tie_bsir_no_prior,12)
expect_identical(names(tie_bsir_no_prior),c("log_posterior","draws","log_proposal","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_length(attributes(tie_bsir_no_prior),12)
expect_identical(names(attributes(tie_bsir_no_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed","where_is_baseline","ncores"))
expect_identical(attr(tie_bsir_no_prior,"approach"),"Bayesian")
expect_silent(print(tie_bsir_no_prior))
expect_silent(summary(tie_bsir_no_prior))
expect_silent(diagnostics(object = tie_bsir_no_prior, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_bsir_no_prior, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_bsir_no_prior,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_error(aic(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(tie_bsir_no_prior)) 


# tests on "WAIC" for "BSIR"                       
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 100L,
                        prior = NULL,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))
tie_reh_bsir_waic <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 100L,
                        prior = NULL,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10)
expect_silent(summary(tie_reh_bsir_waic))

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
                        log = TRUE,
                        seed = set_seed
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
                        log = TRUE,
                        seed = set_seed)

expect_silent(tie_bsir_with_prior)
expect_inherits(tie_bsir_with_prior,"remstimate")  
expect_length(tie_bsir_with_prior,13)
expect_identical(names(tie_bsir_with_prior),c("log_posterior","draws","log_proposal","log_prior","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_length(attributes(tie_bsir_with_prior),13)
expect_identical(names(attributes(tie_bsir_with_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","prior","nsim","seed","where_is_baseline","ncores"))
expect_identical(attr(tie_bsir_with_prior,"approach"),"Bayesian")
expect_silent(print(tie_bsir_with_prior))
expect_silent(summary(tie_bsir_with_prior))
expect_silent(diagnostics(object = tie_bsir_with_prior, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_bsir_with_prior, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_bsir_with_prior,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_error(aic(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(tie_bsir_with_prior)) 

# (4) method  = "HMC"    
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)) 
tie_hmc <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 3L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)
expect_silent(tie_hmc)
expect_inherits(tie_hmc,"remstimate")  
expect_length(tie_hmc,10)
expect_identical(names(tie_hmc),c("draws","log_posterior","coefficients","post.mean","vcov","sd","loglik","df.null","df.model","df.residual"))
expect_length(attributes(tie_hmc),16)
expect_identical(names(attributes(tie_hmc)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed","nchains","burnin","thin","init","where_is_baseline","ncores"))
expect_identical(attr(tie_hmc,"approach"),"Bayesian")
expect_silent(print(tie_hmc))
expect_silent(summary(tie_hmc))
expect_silent(diagnostics(object = tie_hmc, reh = tie_reh, stats = tie_reh_stats))
tie_reh_diagnostics <- diagnostics(object = tie_hmc, reh = tie_reh, stats = tie_reh_stats)
expect_silent(plot(x = tie_hmc,reh = tie_reh, diagnostics = tie_reh_diagnostics))
expect_error(aic(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(tie_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(tie_hmc)) 

# tests on "WAIC" for "HMC"                      
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 100L,
                        burnin = 5L,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))



# testing with omit_dyad
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie", omit_dyad = list(list(time = c(120,148), dyad=data.frame(actor1=c("4",NA),actor2=c(NA,"4"),type=c(NA,NA))))) # removing actor "4" from time=120 to time=148
# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model, method="pt")
# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))

# Risk set "active"

# testing estimation methods with active riskset 
tie_reh <- remify::remify(edgelist = tie_data$edgelist, 
                            model = "tie", 
                            riskset="active") 
# calculating statistics
tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model, method="pt")
attr(tie_reh_stats,"formula") <- NULL

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE"))
tie_mle <- remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = tie_mle, reh = tie_reh, stats = tie_reh_stats))
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
                        prior = NULL,
                        seed = set_seed))     

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
                        log = TRUE,
                        seed = set_seed
                        ))
# (4) method  = "HMC"    
expect_silent(remstimate::remstimate(reh = tie_reh,
                        stats = tie_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)) 



