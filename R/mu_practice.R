#!/usr/bin/env Rscript
#
library(data.table)
library(muStat)

dt <- data.table(
  name = LETTERS[1:6],
  a = c(1, 2, 3, 4, 5, 6),
	b = c(2, 2, 4, 4, 6, 6),
	c = c(9, 8, 7, 7, 8, 9))

# GE = greater than or equal to
ge <- mu.GE(as.matrix(dt[, .(a, b, c)]))

# Logical AND GEs for all variables 
ge_and <- mu.AND(ge)
sums <- mu.Sums(ge_and)
setDT(sums)
sums$name <- dt$name
sums <- setorder(sums, -score)

for (i in 1:nrow(sums)) {
  message(sprintf("RANK: %2d. %s: score = %2d; nAbove = %2d; nBelow = %2d", i, sums$name[i], sums$score[i], sums$nAbove[i], sums$nBelow[i]))
  
}
