
test_that("testing warnings of remstimate", {
  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::inertia()
  
  # calculating statistics
  tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)


  ## seed = c(1234,4321), nchains = 2L, method = "HMC"
  expect_warning(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                         seed = c(1234,4321),
                         nchains = 2L,
                         nsim = 10L,
                         burnin = 5L),
  "`seed` length is greater than 1. Considering only the first element",
  fixed = TRUE)  

  ## seed = c(1234,4321), method = "BSIR"                      
  expect_warning(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "BSIR",
                          nsim = 10L,
                         seed = c(1234,4321)),
  "`seed` length is greater than 1. Considering only the first element",
  fixed = TRUE
  )   

  ## method = "HMC", thin = NULL, nsim >= 100
  expect_warning(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nchains = 2L,
                          nsim = 100L,
                          burnin = 5L,
                          thin = NULL),
  "'thin' parameter undefined. Using thin = 10",
  fixed = TRUE)  

  ## method = "HMC", thin = NULL, nsim >= 100
  expect_warning(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "HMC",
                          nchains = 2L,
                          nsim = 10L,
                          burnin = 5L,
                          thin = NULL),
  "'nsim' is less than 100. No thinning is applied",
  fixed = TRUE)  

})
