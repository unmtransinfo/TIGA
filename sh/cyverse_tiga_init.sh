#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
TIMESTAMP="$(date +'%Y%m%d%H%M%S')"
#
date
###
sudo apt update
sudo apt install -y r-base-dev
sudo R -e 'install.packages(c("readr", "data.table", "igraph", "devtools"))'
#
###
git clone https://github.com/unmtransinfo/TIGA.git
#
###
# SRCDATADIR="/home/kasm-user/data-store/data/swcactiZone/home/jjyang/analyses/GWASCatalog/releases/2026/07/10/"
ln -s $HOME/data-store/data/swcactiZone/home/jjyang/analyses $HOME/data
#
###
mkdir -p $HOME/venv/bioclients
python3 -m venv $HOME/venv/bioclients --clear
source $HOME/venv/bioclients/bin/activate
pip install --upgrade pip
pip install --upgrade bioclients
#
###
date
printf "DONE (${SCRIPTNAME})\n"
###
# Copy VICE log to output folder.
cp $HOME/.vice-init.log $HOME/data-store/data/output/vice-init_${TIMESTAMP}.log
#
