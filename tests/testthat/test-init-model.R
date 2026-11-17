test_that("initial outcome models default to the CATE model for each outcome", {
  b <- simulate_binary(n = 500, seed = 31)
  count <- simulate_count(n = 500, seed = 32)
  surv <- simulate_survival(n = 500, seed = 33)
  cases <- list(
    list(response = "binary", data = b, cate.model = y ~ age,
         init.model = y ~ age, initial = NULL),
    list(response = "count", data = count, cate.model = y ~ age + offset(log(years)),
         init.model = y ~ age + offset(log(years)), initial = "poisson"),
    list(response = "survival", data = surv, cate.model = survival::Surv(y, d) ~ age,
         init.model = survival::Surv(y, d) ~ age, initial = "logistic")
  )
  for (case in cases) {
    args <- list(response = case$response, data = case$data, score.method = "twoReg",
                 cate.model = case$cate.model, ps.model = trt ~ age,
                 prop.cutoff = c(0.5, 1), B = 1, Kfold = 2, seed = 56)
    if (!is.null(case$initial)) args$initial.predictor.method <- case$initial
    if (case$response == "survival") args$tau0 <- median(surv$y) * 1.1
    implicit <- suppressWarnings(do.call(catefit, args))
    explicit <- suppressWarnings(do.call(catefit, c(args, list(init.model = case$init.model))))
    expect_equal(implicit$score.twoReg, explicit$score.twoReg, info = case$response)
    wrong <- if (case$response == "survival") survival::Surv(y + 1, d) ~ age else I(y + 1) ~ age
    expect_error(do.call(catefit, c(args, list(init.model = wrong))), "same outcome", info = case$response)
  }
})

test_that("separate initial designs are used for fitting and validation", {
  for (response in c("binary", "count", "survival")) {
    d <- switch(response, binary = simulate_binary(n = 600, seed = 41),
                count = simulate_count(n = 600, seed = 41),
                survival = simulate_survival(n = 600, seed = 41))
    cm <- switch(response, survival = survival::Surv(y, d) ~ age,
                 count = y ~ age + offset(log(years)), y ~ age)
    im <- switch(response, survival = survival::Surv(y, d) ~ age + female,
                 count = y ~ age + female, y ~ age + female)
    args <- list(response = response, data = d, score.method = "twoReg",
                 cate.model = cm, init.model = im, ps.model = trt ~ age,
                 prop.cutoff = c(0.5, 1), B = 1, Kfold = 2, seed = 23)
    if (response != "binary") args$initial.predictor.method <- switch(response, count = "poisson", survival = "logistic", continuous = "gaussian")
    if (response == "survival") args$tau0 <- median(d$y) * 1.1
    fit <- suppressWarnings(do.call(catefit, args))
    expect_true(all(is.finite(fit$score.twoReg)), info = response)
    same <- args; same$init.model <- cm
    baseline <- suppressWarnings(do.call(catefit, same))
    expect_false(isTRUE(all.equal(fit$score.twoReg, baseline$score.twoReg)), info = response)
    args$cv.n <- 1; args$abc <- FALSE
    if (response != "binary") { args$error.max <- 1; args$max.iter <- 2 }
    cv <- suppressWarnings(do.call(catecv, args))
    expect_true(any(is.finite(cv$ate.twoReg$ate.est.valid.high.cv)), info = response)
  }
})


test_that("continuous preprocessing defaults the initial design and rejects different outcomes", {
  d <- simulate_binary(n = 100, seed = 66)
  args <- list(fun = "catefit", cate.model = y ~ age, init.model = NULL,
               ps.model = trt ~ age, data = d, prop.cutoff = c(0.5, 1),
               ps.method = "glm", score.method = "twoReg")
  out <- do.call(precmed:::data.preproc.mean, args)
  expect_equal(out$x.init, out$x.cate)
  args$init.model <- y ~ age + female
  expect_equal(ncol(do.call(precmed:::data.preproc.mean, args)$x.init), 2)
  args$init.model <- I(y + 1) ~ age
  expect_error(do.call(precmed:::data.preproc.mean, args), "same outcome")
})

test_that("contrast regression accepts separate initial outcome models", {
  for (response in c("count", "survival")) {
    d <- if (response == "count") simulate_count(n = 700, seed = 71) else simulate_survival(n = 700, seed = 71)
    cm <- if (response == "count") y ~ age + offset(log(years)) else survival::Surv(y, d) ~ age
    im <- if (response == "count") y ~ age + female else survival::Surv(y, d) ~ age + female
    args <- list(response = response, data = d, score.method = "contrastReg",
                 cate.model = cm, init.model = im, ps.model = trt ~ age,
                 initial.predictor.method = if (response == "count") "poisson" else "logistic",
                 prop.cutoff = c(0.5, 1), B = 1, Kfold = 2, seed = 72)
    if (response == "survival") args$tau0 <- median(d$y) * 1.1
    fit <- suppressWarnings(do.call(catefit, args))
    expect_true(all(is.finite(fit$score.contrastReg)), info = response)
    args$cv.n <- 1; args$abc <- FALSE; args$error.max <- 1; args$max.iter <- 2
    cv <- suppressWarnings(do.call(catecv, args))
    expect_true(any(is.finite(cv$ate.contrastReg$ate.est.valid.high.cv)), info = response)
  }
})

test_that("count initial offsets must agree with the CATE model", {
  d <- simulate_count(n = 200, seed = 81)
  expect_error(precmed:::initial_outcome_design(d, y ~ age + offset(log(years)),
                                              y ~ age + offset(log(years + 1))), "offset")
  expect_error(precmed:::initial_outcome_design(d, y ~ age, y ~ age + female - 1), "intercept")
})
