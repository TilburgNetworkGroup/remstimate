test_that("testing input arguments and tie-oriented modeling", {

  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::inertia()
  # calculating statistics
  tie_reh_stats <- remstats::remstats(edgelist = tie_reh, tie_effects = tie_model)

  ## input `reh` is a data.frame 
  mle_loc <- remstimate::remstimate(reh = attr(tie_reh,"remulate.reh"),
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1L,
                          model = "tie")

  # testing errors 

  ## input `reh` is null
  expect_error(remstimate::remstimate(reh = NULL,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1L),
  "missing 'reh' argument.",
  fixed = TRUE
  )         

  ## input `reh` is not an `reh` object nor a  `data.frame`
  expect_error(remstimate::remstimate(reh = data.matrix(attr(tie_reh,"remulate.reh ")),
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1L),
  "Input 'reh' must be either a 'reh' object (from 'remify' package) or the event sequence (data.frame).",
  fixed = TRUE
  )       

  ## input `method` is not available or it is mistyped
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "mle",
                          ncores = 1L),
  "The `method` specified is not available or it is mistyped.",
  fixed = TRUE
  )

  ## directed = FALSE and model = "actor"
  data(ao_reh) # loading ao data
  rate_model <- ~ 1 + remstats::indegreeSender() # linear predictor for the rate model
  choice_model <- ~ remstats::inertia() + remstats::reciprocity() # linear predictror for the choice model
  ## calculate statistics
  ao_reh_stats <- remstats::remstats(edgelist = ao_reh, sender_effects = rate_model, 
  receiver_effects = choice_model)
  attr(ao_reh,"directed") <- FALSE
  expect_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          method = "MLE",
                          ncores = 1L),
  "actor-oriented modeling can't operate on undirected networks",                      
  fixed = TRUE
  )
  ## tests on the stats object

  # [ ... here ... ]

  ## input `epochs` is negative
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1L,
                          epochs = -29),
  "'epoch' must be a positive number",
  fixed = TRUE
  )  

  ## input `epochs` is a character
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1L,
                          epochs = "ss"),
  "'epoch' must be a positive number",
  fixed = TRUE
  )  

  ## input `epsilon` must be a positive number
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1L,
                          epsilon = -3.2),
  "'epsilon' must be a positive number",
  fixed = TRUE
  )

  # estimate with NULL epochs and epsilon
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1L,
                          epochs = NULL,
                          epsilon = NULL))


  ## when parallel::detectCores() == 2 and input `ncores` is greater than 1
  if(parallel::detectCores() <= 2L)
  { # run the test only if there are 2 cores in the machine
    expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 10L),
  "'ncores' is recommended to be set at most to 1.",
  fixed = TRUE
  )
  }
  else{
  ## when parallel::detectCores() > 2 and the input `ncores` is greater than the suggested maximum value 
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 2e03L),
  "'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)",
  fixed = TRUE
  )

  
  }

  # estimate with NULL ncores 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = NULL))

  # tests on seed, nchains, nsim, burnin and thin parameters

  ## seed == NULL, nchains = 3L and method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = NULL,
                         nchains = 3L,
                         nsim = 100L,
                         burnin = 5L))
  ## seed == NULL, nchains = NULL and method = "HMC"                      
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = NULL,
                         nchains = NULL,
                         nsim = 100L,
                         burnin = 5L)) 
  ## seed = c(1234,4321), nchains = 2L, method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = 2L,
                         nsim = 100L,
                         burnin = 5L))          

  ## seed = c(1234,4321), nchains = NULL, method = "HMC"   
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = NULL,
                         nsim = 100L,
                         burnin = 5L))      

  ## seed = c(1234,4321), nchains = 5L, method = "HMC"   
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = 5L,
                         nsim = 100L,
                         burnin = 5L),
  "the number of chains (`nchains`) must be equal to the number of seeds (`seed`) supplied",
  fixed = TRUE
  )  

  ## seed = NULL, method = "BSIR" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                          seed = NULL,
                          nsim = 100L)) 

  ## seed = c(1234,4321), method = "BSIR"                      
  expect_warning(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                         seed = c(1234,4321)),
  "`seed` length is greater than 1. Considering only the first element",
  fixed = TRUE
  )   

  ## nsim = NULL, method  = "BSIR" (or "HMC")                    
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                          nsim = NULL,
                          seed = 1234)) 
  ## burnin = NULL, method = "HMC" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          burnin = NULL,
                          seed = 1234)) 

  ## method = "HMC", burnin > nsim      
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          burnin = 1e05L),
  "`burnin` value must be lower than the number of simulations (`nsim`)",
  fixed = TRUE
  )          

  ## burnin = 500L, nsim = 1e03L, thin = NULL, method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          burnin = 500L,
                          nsim = 1e03L,
                          thin = NULL)) 

  ## if((nsim-burnin)<500L), method = "HMC" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 100L,
                          burnin = 50L,
                          thin = NULL))   

  ## if((nsim-burnin)>=500L), method = "HMC" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 600L,
                          burnin = 50L,
                          thin = NULL))   
  ## thin = NULL, method == "HMC" 

  ## thin <= 0, method == "HMC"                               
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 600L,
                          burnin = 50L,
                          thin = 0),
  "`thin` value must be positive. Set thin = 1 if no thinning of chains is required",
  fixed = TRUE
  ) 
})

test_that("testing tie-oriented modeling (methods)", {

  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 

  # calculating statistics
  tie_reh_stats <- remstats::remstats(edgelist = tie_reh, tie_effects = tie_model)

  # (1) method = "MLE"
  tie_mle <- expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1L,
                          method = "MLE"))
  expect_no_error(tie_mle)
  expect_no_error(print(tie_mle))
  expect_no_error(summary(tie_mle))
  expect_no_error(print(summary(tie_mle)))
  expect_no_error(predict(tie_mle))
  expect_no_error(plot(tie_mle))
  expect_no_error(aic(tie_mle))
  expect_no_error(aicc(tie_mle))
  expect_no_error(bic(tie_mle))
  expect_no_error(waic(tie_mle))
  
  # (2) method = "GDADAMAX" 
  tie_gdadamax <- expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1L,
                          method = "GDADAMAX",
                          epochs = 10L))
  expect_no_error(tie_gdadamax)
  expect_no_error(print(tie_gdadamax))
  expect_no_error(summary(tie_gdadamax))
  expect_no_error(print(summary(tie_gdadamax)))
  expect_no_error(predict(tie_gdadamax))
  expect_no_error(plot(tie_gdadamax))
  expect_no_error(aic(tie_gdadamax))
  expect_no_error(aicc(tie_gdadamax))
  expect_no_error(bic(tie_gdadamax))
  expect_no_error(waic(tie_gdadamax))                   

  # (3) method  = "BSIR" 
  
  ## (3.1) with prior = NULL
  tie_bsir_no_prior <- expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1L,
                          method = "BSIR",
                          nsim = 100L,
                          prior = NULL))
  expect_no_error(tie_bsir_no_prior)
  expect_no_error(print(tie_bsir_no_prior))
  expect_no_error(summary(tie_bsir_no_prior))
  expect_no_error(print(summary(tie_bsir_no_prior)))
  expect_no_error(predict(tie_bsir_no_prior))
  expect_no_error(plot(tie_bsir_no_prior))
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

  # (4) method  = "HMC"     
  tie_hmc <- expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1L,
                          method = "HMC",
                          nchains = 2L,
                          nsim = 100L,
                          burnin = 50L))
  expect_no_error(tie_hmc)
  expect_no_error(print(tie_hmc))
  expect_no_error(summary(tie_hmc))
  expect_no_error(print(summary(tie_hmc)))
  expect_no_error(predict(tie_hmc))
  expect_no_error(plot(tie_hmc))
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

})

test_that("testing actor-oriented modeling", {

  # loading data
  data(ao_reh)

  # specifying linear predictor (for rate and choice model)
  rate_model <- ~ 1 + remstats::indegreeSender()
  choice_model <- ~ remstats::inertia() + remstats::reciprocity()

  # calculating statistics
  ao_reh_stats <- remstats::remstats(edgelist = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)

  # (1) method = "MLE"
  ao_mle <- expect_no_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          ncores = 1L,
                          method = "MLE"))
  expect_no_error(ao_mle)
  expect_no_error(print(ao_mle))
  expect_no_error(summary(ao_mle))
  expect_no_error(print(summary(ao_mle)))
  expect_no_error(predict(ao_mle))
  expect_no_error(plot(ao_mle))
  expect_no_error(aic(ao_mle))
  expect_no_error(aicc(ao_mle))
  expect_no_error(bic(ao_mle))
  expect_no_error(waic(ao_mle))

  # (2) method = "GDADAMAX" 
  ao_gdadamax <- expect_no_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          ncores = 1L,
                          method = "GDADAMAX",
                          epochs = 10L))
  expect_no_error(ao_gdadamax)
  expect_no_error(print(ao_gdadamax))
  expect_no_error(summary(ao_gdadamax))
  expect_no_error(print(summary(ao_gdadamax)))
  expect_no_error(predict(ao_gdadamax))
  expect_no_error(plot(ao_gdadamax))
  expect_no_error(aic(ao_gdadamax))
  expect_no_error(aicc(ao_gdadamax))
  expect_no_error(bic(ao_gdadamax))
  expect_no_error(waic(ao_gdadamax)) 

  # (3) method  = "BSIR" 
  
  ## (3.1) with prior = NULL
  ao_bsir_no_prior <- expect_no_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          ncores = 1L,
                          method = "BSIR",
                          nsim = 100L,
                          prior = NULL))
  expect_no_error(ao_bsir_no_prior)
  expect_no_error(print(ao_bsir_no_prior))
  expect_no_error(summary(ao_bsir_no_prior))
  expect_no_error(print(summary(ao_bsir_no_prior)))
  expect_no_error(predict(ao_bsir_no_prior))
  expect_no_error(plot(ao_bsir_no_prior))
  expect_error(aic(ao_bsir_no_prior),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(aicc(ao_bsir_no_prior),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(bic(ao_bsir_no_prior),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(waic(ao_bsir_no_prior),
  "'approach' must be 'Frequentist'",
  fixed = TRUE) 

  ## (3.2) with a specified prior

  # (4) method  = "HMC"     
  ao_hmc <- expect_no_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          ncores = 1L,
                          method = "HMC",
                          nchains = 2L,
                          nsim = 100L,
                          burnin = 50L))
  expect_no_error(ao_hmc)
  expect_no_error(print(ao_hmc))
  expect_no_error(summary(ao_hmc))
  expect_no_error(print(summary(ao_hmc)))
  expect_no_error(predict(ao_hmc))
  expect_no_error(plot(ao_hmc))
  expect_error(aic(ao_hmc),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(aicc(ao_hmc),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(bic(ao_hmc),
  "'approach' must be 'Frequentist'",
  fixed = TRUE)
  expect_error(waic(ao_hmc),
  "'approach' must be 'Frequentist'",
  fixed = TRUE) 


})


test_that("testing warning and errors from methods of a remstimate object", {
  # dummy test here
  expect_equal(2 * 2, 4)
})
