Pareto_filter <- function(dt, col_a, col_b, n) {
  if (nrow(dt) < n) { return(dt[, ok := T]) }
  dt$ok <- F
  n_ok_previous <- 0
  # Include non-dominated solutions up to (n).
  while (sum(dt$ok) < n) {
    if (sum(!dt$ok)==0) { break } #None left -- should not happen.
    i_notok <- which(!dt$ok)
    for (i in i_notok) {
      a <- dt[[col_a]][i]
      b <- dt[[col_b]][i]
      if (sum((!dt$ok) & (dt[[col_a]]>a) & (dt[[col_b]]>b), na.rm=T)==0) {
        dt[i]$ok <- T
        break
      }
    }
    if (sum(dt$ok)==n_ok_previous) { break } #No more -- should not happen.
    n_ok_previous <- sum(dt$ok)
  }
  return(dt)
}
