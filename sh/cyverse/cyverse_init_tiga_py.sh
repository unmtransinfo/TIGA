#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
#
date
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
