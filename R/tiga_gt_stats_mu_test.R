library(readr)
library(data.table)
gt_stats <- read_delim("data/20201216/gt_stats_mu.tsv.gz", "\t")
setDT(gt_stats)

message(sprintf("nAbove: max: %.1f; min: %.1f\n\tquantiles:", min(gt_stats$nAbove), max(gt_stats$nAbove)))
print(quantile(gt_stats$nAbove))
message(sprintf("nBelow: max: %.1f; min: %.1f\n\tquantiles:", min(gt_stats$nBelow), max(gt_stats$nBelow)))
print(quantile(gt_stats$nBelow))
message(sprintf("muScore: max: %.1f; min: %.1f\n\tquantiles:", min(gt_stats$muScore), max(gt_stats$muScore)))
print(quantile(gt_stats$muScore))

setorder(gt_stats, -muScore)
gt_stats[, muRank := 1:nrow(gt_stats)]

