#!/bin/bash
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
#
VEGAS="${cwd}/perl/vegas2v2.pl"
#
cd $DATADIR
#
$VEGAS -G \
	-snpandp /home/data/VEGAS/data/VEGAS2v2example/example.txt \
	-custom /home/data/VEGAS/data/VEGAS2v2example/example \
	-glist /home/data/VEGAS/data/VEGAS2v2example/example.glist \
	-out vegas_example
#
