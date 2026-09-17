#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
TIMESTAMP="$(date +'%Y%m%d%H%M%S')"
#
date
###
# SRCDATADIR="/home/kasm-user/data-store/data/swcactiZone/home/jjyang/analyses/GWASCatalog/releases/2026/07/10/"
ln -s $HOME/data-store/data/swcactiZone/home/jjyang/analyses $HOME/data
###
$HOME/data/cyverse_init_tiga_r.sh >& $HOME/cyverse_init_tiga_r_${TIMESTAMP}.log
cp $HOME/cyverse_init_tiga_r_${TIMESTAMP}.log $HOME/data-store/data/output/
###
#
###
$HOME/data/cyverse_init_tiga_py.sh >& $HOME/cyverse_init_tiga_py_${TIMESTAMP}.log
cp $HOME/cyverse_init_tiga_py_${TIMESTAMP}.log $HOME/data-store/data/output/
###
git clone https://github.com/unmtransinfo/TIGA.git
#
cd TIGA
./sh/Go_TIGA_Workflow.sh >& $HOME/Go_TIGA_Workflow_${TIMESTAMP}.log
cp $HOME/Go_TIGA_Workflow_${TIMESTAMP}.log $HOME/data-store/data/output/
###
date
printf "DONE (${SCRIPTNAME})\n"
###
# Copy VICE log to output folder.
cp $HOME/.vice-init.log $HOME/data-store/data/output/vice-init_${TIMESTAMP}.log
#
