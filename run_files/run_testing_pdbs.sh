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
#echo "User: $(whoami)" # Fails 
echo "Working Directory: $(pwd)"
echo "============================="
echo ""

echo "=== Initial Files in Directory ==="
ls -la

echo "============================="
echo "PROGRAM RUNNING WITH " ${PDB}
echo "============================="

echo "=== Initial Files in SuperWater Dir ==="
ls /SuperWater -ltr

pip install MDAnalysis

echo "starting download"
python download_list.py test_res15.txt

echo "# downloaded this many PDB files:" 
ls *pdb | wc
echo "# which takes up this much space"
du -sh

pwd
echo "# running superwaters"
superwater-predict --config ./super_train_config.yaml --verbose
echo "# finished; taring ..."
tar -czf supa_outs_${PDB}_${CLUSTER}.tar.gz ./supa_outs
echo "# finished, removing .pdb files"

# bypass argument length and will remove child dir
find . -name "*.pdb" -print0 | xargs -0 rm
# rm sub dir also
rm -r ./supa_outs/

echo "KYRIE ELESON"
