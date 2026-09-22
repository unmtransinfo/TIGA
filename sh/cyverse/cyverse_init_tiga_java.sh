#!/bin/bash
#
###
SCRIPTNAME="$(basename $0)"
#
date
###
sudo apt update
#
sudo apt install -y openjdk-21-jdk
sudo apt install -y maven
#
mkdir -p ~/app/lib
#
cd ~/app
git clone https://github.com/IUIDSL/iu_idsl_jena.git
cd iu_idsl_jena
#
mvn compile
mvn install
cp target/*.jar ~/app/lib
###
date
printf "DONE (${SCRIPTNAME})\n"
#
