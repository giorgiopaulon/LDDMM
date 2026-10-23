# The label update is checked by enumeration: starting from exact posterior
# draws, one update must leave the distribution unchanged.

test_that("the label update leaves the exact posterior invariant", {
  skip_if_not_installed("statmod")
  for (flexible in c(FALSE, TRUE)){
    p = toy_label_problem(flexible)
    probs = toy_exact(p)
    for (informed in c(TRUE, FALSE)){
      set.seed(2)
      counts = toy_one_step(p, probs, 20000, informed)
      expect_gt(chisq_pvalue(counts, probs), 1e-3)
    }
  }
})

test_that("approximate kernels fail the same check (the test has power)", {
  skip_if_not_installed("statmod")
  p = toy_label_problem(TRUE)
  probs = toy_exact(p)
  set.seed(2)
  expect_lt(chisq_pvalue(toy_one_step(p, probs, 20000, TRUE, drop_norm = TRUE), probs), 1e-6)
  expect_lt(chisq_pvalue(toy_one_step(p, probs, 20000, TRUE, simple_prior = TRUE), probs), 1e-6)
})
