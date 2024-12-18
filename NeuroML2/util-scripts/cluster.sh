#!/bin/bash

# Copyright 2024 Ankur Sinha
# Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
# File : cluster.sh
#
# Utility script to do things with the cluster.

CLUSTER_SSH="ucl-myriad"
CLUSTER_SIM_PATH="./Scratch/HL23/"
SIM_FOLDER=""


send_sim_to_cluster () {
    echo "> Pushing simulation files to cluster"
    rsync -avPh ${SIM_FOLDER} ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}
}

get_logs_data () {
    echo "> Getting logs and data files"
    rsync -avPh ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}/${SIM_FOLDER}/*.sh.* ${SIM_FOLDER}/
    rsync -avPh ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}/${SIM_FOLDER}/*.dat ${SIM_FOLDER}/
}

usage () {
    echo "$0: utility script to interact with cluster"
    echo
    echo "-s <simulation folder>"
    echo "   put simulation on cluster"
    echo "-g <simulation folder>"
    echo "   get simulation logs and data from cluster"
    echo "-h: print help and exit"
}

while getopts "s:g:h" OPTION
do
    case $OPTION in
        s)
            SIM_FOLDER="$OPTARG"
            send_sim_to_cluster
            exit 0
            ;;
        g)
            SIM_FOLDER="$OPTARG"
            get_logs_data
            exit 0
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 0
            ;;
    esac
done
