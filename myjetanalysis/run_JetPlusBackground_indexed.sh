#!/bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME} 
 
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/verkest/install
 
# print the environment - needed for debugging
printenv
echo 1: $1

root.exe -q -b macro/Fun4All_JetPlusBackground_forCondor.C\($1,$2,$3,$4\)

echo all done
