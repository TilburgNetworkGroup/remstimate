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
                          ncores = 1,
                          model = "tie")

  # testing errors 

  ## input `reh` is null
  expect_error(remstimate::remstimate(reh = NULL,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1),
  "missing 'reh' argument.",
  fixed = TRUE
  )         

  ## input `reh` is not an `reh` object nor a  `data.frame`
  expect_error(remstimate::remstimate(reh = data.matrix(attr(tie_reh,"remulate.reh ")),
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 1),
  "Input 'reh' must be either a 'reh' object (from 'remify' package) or the event sequence (data.frame).",
  fixed = TRUE
  )       

  ## input `method` is not available or it is mistyped
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "mle",
                          ncores = 1),
  "The `method` specified is not available or it is mistyped.",
  fixed = TRUE
  )  

  ## tests on the stats object

  # [ ... here ... ]

  ## input `epochs` is negative
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1,
                          epochs = -29),
  "'epoch' must be a positive number",
  fixed = TRUE
  )  

  ## input `epochs` is a character
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1,
                          epochs = "ss"),
  "'epoch' must be a positive number",
  fixed = TRUE
  )  

  ## input `epsilon` must be a positive number
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1,
                          epsilon = -3.2),
  "'epsilon' must be a positive number",
  fixed = TRUE
  )

  # estimate with NULL epochs and epsilon
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "GDADAMAX",
                          ncores = 1,
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
                         nsim = 100,
                         burnin = 5))
  ## seed == NULL, nchains = NULL and method = "HMC"                      
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = NULL,
                         nchains = NULL,
                         nsim = 100,
                         burnin = 5)) 
  ## seed = c(1234,4321), nchains = 2L, method = "HMC"
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = 2L,
                         nsim = 100,
                         burnin = 5))          

  ## seed = c(1234,4321), nchains = NULL, method = "HMC"   
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = NULL,
                         nsim = 100,
                         burnin = 5))      

  ## seed = c(1234,4321), nchains = 5L, method = "HMC"   
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = 5L,
                         nsim = 100,
                         burnin = 5),
  "the number of chains (`nchains`) must be equal to the number of seeds (`seed`) supplied",
  fixed = TRUE
  )  

  ## seed = NULL, method = "BSIR" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                          seed = NULL,
                          nsim = 100)) 

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
                          method = "HMC")) 

  ## if((nsim-burnin)<500L), method = "HMC" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 100,
                          burnin = 50))   

  ## if((nsim-burnin)>=500L), method = "HMC" 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nsim = 600,
                          burnin = 50))                 

})

test_that("testing tie-oriented modeling (methods)", {

  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 

  # calculating statistics
  tie_reh_stats <- remstats::remstats(edgelist = tie_reh, tie_effects = tie_model)

  # (1) method = "MLE"
  tie_mle <- remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1,
                          method = "MLE") 
  # these two tests are risky because estimates can slightly change across operating systems                          
  # expect_snapshot(tie_mle)
  # expect_snapshot(summary(tie_mle))
  
  # dummy test here
  expect_equal(2 * 2, 4)
})

test_that("testing actor-oriented modeling", {

  # dummy test here
  expect_equal(2 + 2, 4)
})
