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
SCRIPTNAME_R="cyverse_init_tiga_r.sh"
LOG="${SCRIPTNAME_R}-${TIMESTAMP}.log"
setsid nohup $WORKDIR/$SCRIPTNAME_R >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPTNAME_R started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
sleep 3
###
SCRIPTNAME_PY="cyverse_init_tiga_py.sh"
LOG="${SCRIPTNAME_PY}-${TIMESTAMP}.log"
setsid nohup $WORKDIR/$SCRIPTNAME_PY >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPTNAME_PY started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
sleep 3
###
SCRIPTNAME_JAVA="cyverse_init_tiga_java.sh"
LOG="${SCRIPTNAME_JAVA}-${TIMESTAMP}.log"
setsid nohup $WORKDIR/$SCRIPTNAME_JAVA >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPTNAME_JAVA started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
sleep 3
###
cd $HOME
git clone https://github.com/unmtransinfo/TIGA.git
###
# TIGA output files to iRODS output dir:
ln -s $OUTDIR/tiga_data $HOME/TIGA/data
#
cp $HOME/data/.tcrd.yaml $HOME/
#
###
# Wait until background init tasks done.
while [ 1 ]; do
	if [ -e "$HOME/${SCRIPTNAME_PY}_DONE.txt" \
		-a -e "$HOME/${SCRIPTNAME_JAVA}_DONE.txt" \
		-a -e "$HOME/${SCRIPTNAME_R}_DONE.txt" ]; then
		break
	else
		sleep 60
	fi
done
#
cd $HOME/TIGA
SCRIPTNAME_TIGA="Go_TIGA_Workflow.sh"
LOG="${SCRIPTNAME_TIGA}-${TIMESTAMP}.log"
./sh/$SCRIPTNAME_TIGA >$OUTDIR/$LOG 2>&1 & disown
echo "$SCRIPTNAME_TIGA started in background (pid $!); tail -f $OUTDIR/$LOG to watch progress"
###
date
printf "DONE (${SCRIPTNAME})\n"
###
# Copy VICE log to output folder.
cp $HOME/.vice-init.log $HOME/data-store/data/output/vice-init_${TIMESTAMP}.log
#
