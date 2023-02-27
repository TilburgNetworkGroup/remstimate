test_that("testing input arguments", {

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
  if(parallel::detectCores() ==2)
  { # run the test only if there are 2 cores in the machine
    expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 10),
  "'ncores' is recommended to be set at most to 1.",
  fixed = TRUE
  )
  }

  ## when parallel::detectCores() > 2 and the input `ncores` is greater than the suggested maximum value 
  expect_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 20),
  paste("'ncores' is recommended to be set at most to",floor(parallel::detectCores()-2),".",sep=" "),
  fixed = TRUE
  )

  # estimate with NULL ncores 
  expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = NULL))

 

})

test_that("testing tie-oriented modeling", {

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
