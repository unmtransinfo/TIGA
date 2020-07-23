#!/usr/bin/env Rscript
###
# mu.GE: greater than or equal to
# mu.AND: Logical AND GEs for all variables 
#
library(data.table)
library(muStat)

###
# Toy example
dt <- data.table(
  name = LETTERS[1:6],
	a = c(1, 2, 3, 4, 5, 6),
	b = c(2, 2, 4, 4, 6, 6),
	c = c(9, 8, 7, 7, 8, 9))

ge <- mu.GE(as.matrix(dt[, .(a, b, c)]))
ge_and <- mu.AND(ge)
sums <- mu.Sums(ge_and)
setDT(sums)
sums$name <- dt$name
sums <- setorder(sums, -score)

for (i in 1:nrow(sums)) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
}

###
# Bigger example (100x7)
n <- 300
dt <- data.table(
	names = sprintf("case_%03d", 1:n),
	a = as.integer(rnorm(n, sd=5)),
	b = as.integer(rnorm(n, sd=5)),
	c = as.integer(rnorm(n, sd=5)),
	d = as.integer(rnorm(n, sd=5)),
	e = as.integer(rnorm(n, sd=5)),
	f = as.integer(rnorm(n, sd=5)),
	g = as.integer(rnorm(n, sd=5)))

ge <- mu.GE(as.matrix(dt[, .(a, b, c, d, e, f, g)]))
ge_and <- mu.AND(ge)
sums <- mu.Sums(ge_and)
setDT(sums)
sums$name <- dt$name
sums <- setorder(sums, -score, na.last=T)

for (i in c(1:10, (nrow(dt)-10):nrow(dt))) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
}

# Handling NAs.
#dt[runif(n)<.1, a := NA]
#dt[runif(n)<.1, f := NA]
