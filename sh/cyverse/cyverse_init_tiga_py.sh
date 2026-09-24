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
touch $HOME/${SCRIPTNAME}_DONE.txt
date 2>&1 >>$HOME/${SCRIPTNAME}_DONE.txt
printf "DONE (${SCRIPTNAME})\n" 2>&1 >>$HOME/${SCRIPTNAME}_DONE.txt
#
date
printf "DONE (${SCRIPTNAME})\n"
###
