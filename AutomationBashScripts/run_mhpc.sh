#!/bin/bash


echo "------------Running MHPC------------"

MHPC_env='/home/lab/Desktop/Nganga/MPCStack/build/'
cd "$MHPC_env"
./MHPC/mhpc_run &
RUNNING_PID=$!

#Kill the process
cleanup() {
    echo "Killing the MHPC bash..."
    kill $RUNNING_PID
    wait $RUNNING_PID
    exit 0
}

#signals to kill the procs 
trap cleanup SIGTERM SIGINT
#wait for resources release
wait
