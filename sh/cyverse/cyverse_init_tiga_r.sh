#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
#
date
###
sudo apt update
sudo apt install -y r-base-dev
sudo R -e 'install.packages(c("readr", "data.table", "igraph", "devtools"))'
#
###
date
printf "DONE (${SCRIPTNAME})\n"
#
