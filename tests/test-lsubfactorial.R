stopifnot(all.equal(
  exp(stat.extend:::lsubfactorial(1:7)),
  c(0, 1, 2, 9, 44, 265, 1854) # verified in wolfram alpha
))

stopifnot(all.equal(
  exp(stat.extend:::lsubfactorial(c(0, NA, 2, -1, 3))),
  c(1, NA, 1, NA, 2)
))