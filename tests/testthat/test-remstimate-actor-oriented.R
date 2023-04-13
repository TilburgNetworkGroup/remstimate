test_that("testing actor-oriented modeling", {

  # loading data
  data(ao_reh)

  # specifying linear predictor (for rate and choice model)
  rate_model <- ~ 1 + remstats::indegreeSender()
  choice_model <- ~ remstats::inertia() + remstats::reciprocity()

  # calculating statistics
  ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)

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
                          epochs = 10,
                          epsilon = NULL))
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
                          nsim = 10L,
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
                          nchains = 1L,
                          nsim = 10L,
                          burnin = 5L))
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