test_that("testing input arguments", {

  # specifying linear predictor
  tie_model <- ~ 1 + remstats::inertia()
  # calculating statistics
  tie_reh_stats <- remstats::remstats(edgelist = tie_reh, tie_effects = tie_model)

  ## input `reh` is a data.frame 
  mle_loc <- remstimate::remstimate(reh = attr(tie_reh,"remulate.reh"),
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 7,
                          model = "tie")

  # testing errors 

  ## input `reh` is null
  expect_error(remstimate::remstimate(reh = NULL,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 7),
  "missing 'reh' argument.",
  fixed = TRUE
  )         

  ## input `reh` is not an `reh` object nor a  `data.frame`
             


})

test_that("testing tie-oriented modeling", {
  # specifying linear predictor
  tie_model <- ~ 1 + remstats::indegreeSender()+remstats::inertia()+remstats::reciprocity() 
  # calculating statistics
  tie_reh_stats <- remstats::remstats(edgelist = tie_reh, tie_effects = tie_model)

  # (1) method = "MLE"
  tie_mle <- remstimate::remstimate(reh = tie_reh,
                          stats = tie_reh_stats,
                          method = "MLE",
                          ncores = 7) 
  # these two tests are risky because estimates can slightly change across operating systems                          
  expect_snapshot(tie_mle)
  expect_snapshot(summary(tie_mle))
})

test_that("testing actor-oriented modeling", {

  # dummy test here
  expect_equal(2 + 2, 4)
})
