#!/bin/bash

# Copyright 2024 Ankur Sinha
# Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
# File : cluster.sh
#
# Utility script to do things with the cluster.
# Note that this script can only be used from within UCL and assumes that you
# have your ssh keys set up for passwordless ssh access.

CLUSTER_SSH="ucl-myriad"
CLUSTER_SIM_PATH="./Scratch/HL23/"
CLUSTER_USERNAME=""
SIM_FOLDER=""


send_sim_to_cluster () {
    echo "> Pushing simulation files to cluster"
    echo "> Command: rsync -avPh ${SIM_FOLDER} ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}"
    rsync -avPh ${SIM_FOLDER} ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}
}

get_logs_data () {
    echo "> Getting logs and data files"
    echo "> Command: rsync -avPh ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}/${SIM_FOLDER}/${SIM_FOLDER}_generated/*.{sh.*,dat} ${SIM_FOLDER}/${SIM_FOLDER}_generated"
    echo
    rsync -avPh ${CLUSTER_SSH}:./${CLUSTER_SIM_PATH}/${SIM_FOLDER}/${SIM_FOLDER}_generated/*.{sh.*,dat} ${SIM_FOLDER}/${SIM_FOLDER}_generated
}

queue_simulation () {
    echo "> Queing simulation"
    echo "> Command: ssh -v ${CLUSTER_SSH} \"{ cd ~/${CLUSTER_SIM_PATH}/${SIM_FOLDER}/${SIM_FOLDER}_generated/ && qsub LEMS*.xml.sh ; }\""
    echo
    ssh -v ${CLUSTER_SSH} "{ cd ~/${CLUSTER_SIM_PATH}/${SIM_FOLDER}/${SIM_FOLDER}_generated/ && qsub LEMS*.xml.sh ; }"
}

get_queue_list () {
    echo "> Queue list:"
    echo "> Command: ssh -v ${CLUSTER_SSH} \"{ qstat -l ${CLUSTER_USERNAME} ; }\""
    echo
    ssh ${CLUSTER_SSH} "{ qstat -l ${CLUSTER_USERNAME} ; }"

}

usage () {
    echo "$0: utility script to interact with cluster"
    echo
    echo "The <simulation folder> is the main folder that will contain another <.._generated> folder in it"
    echo
    echo "-s <simulation folder>"
    echo "   put simulation on cluster"
    echo "-q <simulation folder>"
    echo "   queue simulation"
    echo "-a <simulation folder>"
    echo "   put simulation on cluster and queue it"
    echo "-g <simulation folder>"
    echo "   get simulation logs and data from cluster"
    echo "-l <cluster username>"
    echo "   get cluster queue list"
    echo "-h: print help and exit"
}

while getopts "s:g:q:a:l:h" OPTION
do
    case $OPTION in
        s)
            SIM_FOLDER="${OPTARG%/}"
            send_sim_to_cluster
            exit 0
            ;;
        q)
            SIM_FOLDER="${OPTARG%/}"
            queue_simulation
            exit 0
            ;;
        a)
            SIM_FOLDER="${OPTARG%/}"
            send_sim_to_cluster
            queue_simulation
            exit 0
            ;;
        g)
            SIM_FOLDER="${OPTARG%/}"
            get_logs_data
            exit 0
            ;;
        l)
            CLUSTER_USERNAME="${OPTARG}"
            get_queue_list
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
