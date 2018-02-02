#!/bin/bash

tasks=""

if [ -z ${1} ]; then
    tasks=`ls crab_* | grep crab_ | sed s/:/\ /g`
else
    if [ -f ${1} ]; then
        tasks=`cat ${1}`
    else
        tasks=`ls ${1} | grep crab_ | sed s/:/\ /g`
    fi
fi

for task in $tasks
do
    crab resubmit $task
done
