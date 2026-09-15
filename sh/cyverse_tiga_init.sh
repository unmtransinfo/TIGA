#!/bin/bash
#
###
mkdir -p $HOME/venv/bioclients
python3 -m venv --clear $HOME/venv/bioclients 
source $HOME/venv/bioclients/bin/activate
pip install --upgrade pip
pip install --upgrade bioclients
#
###
sudo apt update
sudo apt install -y r-base-dev
sudo R -e 'install.packages(c("readr", "data.table", "igraph", "devtools"))
#
###
git clone https://github.com/unmtransinfo/TIGA.git
#
###
# SRCDATADIR="/home/kasm-user/data-store/data/swcactiZone/home/jjyang/analyses/GWASCatalog/releases/2026/07/10/"
#
###
ln -s /home/kasm-user/data-store/data/swcactiZone/home/jjyang/analyses $HOME/data
#
