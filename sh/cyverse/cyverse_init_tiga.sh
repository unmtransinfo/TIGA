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
#
WORKDIR="$HOME/data"
OUTDIR="$HOME/data-store/data/output"
###
SCRIPT="cyverse_init_tiga_r.sh"
LOG="${SCRIPT}-${TIMESTAMP}.log"
setsid nohup $WORKDIR/$SCRIPT >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPT started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
sleep 3
###
SCRIPT="cyverse_init_tiga_py.sh"
LOG="${SCRIPT}-${TIMESTAMP}.log"
setsid nohup $WORKDIR/$SCRIPT >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPT started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
sleep 3
###
git clone https://github.com/unmtransinfo/TIGA.git
#
#cd TIGA
#./sh/Go_TIGA_Workflow.sh >& $OUTDIR/Go_TIGA_Workflow_${TIMESTAMP}.log
###
date
printf "DONE (${SCRIPTNAME})\n"
###
# Copy VICE log to output folder.
cp $HOME/.vice-init.log $HOME/data-store/data/output/vice-init_${TIMESTAMP}.log
#
