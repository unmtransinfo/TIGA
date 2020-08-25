#!/usr/bin/env Rscript
###
# mu.GE: greater than or equal to
# mu.AND: Logical AND GEs for all variables 
###
# If a case has nAbove=0 and nBelow=0 then weight=0 and mu.Sums() returns score = NA.
# This is arbitrary so we assign score = 0 - 0 = 0.
# ?mu.Sums:
# score  = (nB-nA) * ifelse(weight==0, NA, 1)
###
#
library(data.table)
library(muStat)


Calculate_mu_scores <- function(dt) {
  #ge <- mu.GE(as.matrix(dt[, .(a, b, c, d, e, f, g)]))
  ge <- mu.GE(as.matrix(dt[, .SD, .SDcols = !c("name")]))
  ge_and <- mu.AND(ge)
  sums <- mu.Sums(ge_and)
  setDT(sums)
  sums$name <- dt$name
  sums[is.na(score) & weight==0, score := 0] # Override NA scores
  sums[order(-score, -weight), rank := 1:nrow(sums)]
  return(sums)
}

###
# Toy example
dt1 <- data.table(
  name = LETTERS[1:6],
	a = c(1, 2, 3, 4, 5, 6),
	b = c(2, 2, 4, 4, 6, 6),
	c = c(9, 8, 7, 7, 8, 9))

sums <- Calculate_mu_scores(dt1)
sums <- setorder(sums, rank)
for (i in 1:nrow(sums)) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
}

###
# Bigger example (100x7)
n <- 100
dt2 <- data.table(
	name = sprintf("case_%03d", 1:n),
	a = as.integer(rnorm(n, sd=5)),
	b = as.integer(rnorm(n, sd=5)),
	c = as.integer(rnorm(n, sd=5)),
	d = as.integer(rnorm(n, sd=5)),
	e = as.integer(rnorm(n, sd=5)),
	f = as.integer(rnorm(n, sd=5)),
	g = as.integer(rnorm(n, sd=5)))

sums <- Calculate_mu_scores(dt2)
sums <- setorder(sums, rank)
for (i in c(1:5, (nrow(sums)-5):nrow(sums))) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
}

# Check handling NAs.
dt3 <- dt2
dt3[runif(n)<.1, a := NA]
dt3[runif(n)<.1, f := NA]

sums <- Calculate_mu_scores(dt3)
sums <- setorder(sums, rank)
for (i in c(1:5, (nrow(sums)-5):nrow(sums))) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
}