#!/usr/bin/env bash

#set -e

mkdir -p output/condor_logs
mkdir -p results

echo "=== HTCondor Job Information ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "System: $(uname -spo)"
echo "_CONDOR_JOB_IWD: $_CONDOR_JOB_IWD"
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"
echo "RunningOn: $RUNNINGON"
echo "User: $(whoami)" # Fails 
#echo "Working Directory: $(pwd)"
echo "============================="
echo ""


echo "============================="
echo "PROGRAM RUNNING WITH " ${PDB} ${RHO}
echo "============================="
python3 run_rism_master_cluster.py ${PDB} ${RHO}
echo "removing files you dont want"

rm *ts4s
rm *gro
rm *log
rm guv*
rm ex*
rm *top
rm *solute
rm *pdb
rm *xz
echo "KYRIE ELESON"
