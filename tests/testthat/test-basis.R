test_that("quadratic B-splines form a partition of unity with weights 1/2 at the knots", {
  knots = 3:9
  x = seq(3, 9, by = 0.05)
  expect_equal(rowSums(B_basis(x, knots)), rep(1, length(x)))
  Bk = B_basis(knots, knots)
  expect_equal(dim(Bk), c(length(knots), length(knots) + 1))
  for (t in seq_along(knots)){
    expect_equal(Bk[t, ], replace(numeric(length(knots) + 1), c(t, t + 1), 0.5))
  }
})

test_that("P_smooth1 is the first-difference penalty", {
  beta = stats::rnorm(6)
  expect_equal(drop(t(beta) %*% P_smooth1(6) %*% beta), sum(diff(beta)^2))
})

test_that("Hamming balls list every vector within the radius exactly once", {
  x = c(1, 3, 2)
  for (r in 1:2){
    ball = H_ball(x, S = 4, r = r)
    expect_equal(ncol(ball), sum(choose(3, 0:r) * 3^(0:r)))
    expect_equal(ball[, 1], x)
    expect_equal(anyDuplicated(t(ball)), 0)
    expect_true(all(colSums(ball != x) <= r))
  }
})
