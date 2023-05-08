context("Families in combTMB")

test_that("Families return a name to list with the correct names", {
  .names <- c("family", "link", "linkfun", "linkinv", "mu.eta")
  expect_true(all(.names %in% names(poigamma(link = "log"))))
  expect_true(all(.names %in% names(betabernoulli_al(link = "probit"))))
})


test_that("The combTMB families work with appropriate links", {
  expect_identical(class(poigamma(link = "log")), "family")
  expect_identical(class(poigamma(link = log)), "family")
  expect_error(class(poigamma(link = "orange")))
  expect_error(class(poigamma(link = orange)))

})
