# Some simple checks on matching

library(stat.extend)

x <- rmatching(1000, 5)

stopifnot(!4 %in% x)


fx <- dmatching(0:5, 5)

stopifnot(fx[4 + 1] == 0)

stopifnot(all.equal(sum(fx), 1))

stopifnot(all.equal(cumsum(fx), pmatching(0:5, 5)))

stopifnot(!is.unsorted(pmatching(qmatching(0:100/100, 5), 5)))


stopifnot(all.equal(
  moments.matching(5),
  structure(list(mean = 1, var = 1, skew = 1, kurt = 4, excess.kurt = 1), row.names = c(NA, -1L), class = "data.frame")
))