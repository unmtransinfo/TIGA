#!/usr/bin/env Rscript
###
# See .Rprofile for port (3838).
###
require(shiny, quietly = T)
#
# /srv/shiny-server/tiga/
#
#	port = getOption("shiny.port"),
runApp(appDir = paste0(Sys.getenv()["HOME"],'/projects/idg/gwas/R/tiga'),
	port = 9999,
	display.mode = "auto", launch.browser = T)
