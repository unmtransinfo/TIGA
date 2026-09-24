#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
#
date
###
sudo apt update
#
# Dependency for shiny:
sudo apt install -y libuv1-dev
# Dependency for plotly:
sudo apt install -y libssl-dev
sudo apt install -y libcurl4-openssl-dev
#
sudo apt install -y r-base-dev
sudo R -e 'install.packages(c("readr", "data.table", "igraph", "shiny", "DT", "shinyBS", "tableHTML", "plotly"))'
#
###
touch $HOME/${SCRIPTNAME}_DONE.txt
date 2>&1 >>$HOME/${SCRIPTNAME}_DONE.txt
printf "DONE (${SCRIPTNAME})\n" 2>&1 >>$HOME/${SCRIPTNAME}_DONE.txt
#
###
date
printf "DONE (${SCRIPTNAME})\n"
#
