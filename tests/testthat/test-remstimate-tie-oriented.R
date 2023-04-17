test_that("testing tie-oriented modeling (methods)", {

  # loading data
  data(tie_reh)

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 

  # calculating statistics
  tie_reh_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

  # (1) method = "MLE"
  tie_mle <- expect_no_error(remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          ncores = 1L,
                          method = "MLE"))
  expect_no_error(tie_mle)
  expect_no_error(print(tie_mle))
  expect_no_error(summary(tie_mle))
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
                          nsim = 10L,
                          prior = NULL))
  expect_no_error(tie_bsir_no_prior)
  expect_no_error(print(tie_bsir_no_prior))
  expect_no_error(summary(tie_bsir_no_prior))
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
                          nchains = 1L,
                          nsim = 10L,
                          burnin = 5L))
  expect_no_error(tie_hmc)
  expect_no_error(print(tie_hmc))
  expect_no_error(summary(tie_hmc))
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
