test_that("testing errors of remstimate", {

  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::inertia()
  
  # calculating statistics
  tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

  ## input `reh` is a data.frame 
  mle_loc <- remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1L,
                          model = "tie")

  # testing errors 

  ## input `reh` is not a `remify` object
  expect_error(remstimate::remstimate(reh = data.matrix(attr(tie_reh,"remulate.reh ")),
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1L),
  "'reh' must be a 'remify' object (see ?remify::remify).",
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
  ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)
  attr(ao_reh,"directed") <- FALSE
  expect_error(remstimate::remstimate(reh = ao_reh,
                          stats = ao_reh_stats,
                          method = "MLE",
                          ncores = 1L),
  "actor-oriented modeling can't operate on undirected networks",                      
  fixed = TRUE
  )
  ## tests on the stats object

  # [ ... HERE ... ]

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

  # estimate with epochs = NULL and epsilon = NULL
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1L,
                          epochs = NULL,
                          epsilon = NULL))


  ## when parallel::detectCores() == 2 and input `ncores` is greater than 1
  if(parallel::detectCores() <= 2L)
  { 
    # run the test only if there are 2 cores in the machine
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

  # tests on seed, nchains, nsim, burnin, thin, L and espilon parameters

  ## seed == NULL, nchains = 3L and method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = NULL,
                         nchains = 1L,
                         nsim = 10L,
                         burnin = 5L))
  ## seed == NULL, nchains = NULL and method = "HMC"                      
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = NULL,
                         nchains = NULL,
                         nsim = 10L,
                         burnin = 5L))         

  ## seed = 1234, nchains = NULL, method = "HMC"   
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = 1234,
                         nchains = NULL,
                         nsim = 10L,
                         burnin = 5L))

  ## seed = NULL, method = "BSIR" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                          seed = NULL,
                          nsim = 10L))

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
                          nsim = 10L,
                          thin = 5L,
                          burnin = NULL,
                          seed = 1234))          

  ##  thin = NULL, method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 10L,
                          burnin = 5L,
                          thin = NULL)) 

  ## thin <= 0, method == "HMC"                               
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 10L,
                          burnin = 5L,
                          thin = 0),
  "`thin` value must be positive. Set thin = 1 if no thinning of chains is required",
  fixed = TRUE
  ) 
})