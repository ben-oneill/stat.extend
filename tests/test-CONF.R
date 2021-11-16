
stopifnot(all.equal(
  stat.extend::CONF.mean(.05, mtcars$mpg, N=1000),
  structure(list(list(l = 17.9527283309983, r = 22.2285216690017, 
                      lc = TRUE, rc = TRUE)), class = c("ci", "interval"), domain = "R", method = NA_character_, data = "Interval uses 32 data points from data mtcars$mpg with sample variance = 36.3241 and assumed kurtosis = 3.0000", confidence = 0.95, parameter = "mean for population of size 1000")))

stopifnot(all.equal(
  stat.extend::CONF.var(.05, mtcars$mpg, N=1000),
  structure(list(list(l = 22.0730002460278, r = 57.9538604068809, 
                      lc = TRUE, rc = TRUE)), class = c("ci", "interval"), domain = "R", method = "Computed using nlm optimisation with 7 iterations (code = 1)", data = "Interval uses 32 data points from data mtcars$mpg with sample variance = 36.3241 and sample kurtosis = 2.7995", confidence = 0.95, parameter = "variance for population of size 1000")
))

stopifnot(all.equal(
  stat.extend::CONF.prop(.05, mtcars$mpg > 20, N=1000),
  structure(list(list(l = 0.283745994666663, r = 0.604278989534686, 
                      lc = TRUE, rc = TRUE)), class = c("ci", "interval"), domain = "R", data = "Interval uses 32 binary data points from data mtcars$mpg > 20 with sample proportion = 0.4375", method = NA_character_, confidence = 0.95, parameter = "proportion for population of size 1000")
))
