#!/bin/bash

popmodel=$1
numsamps=$2
priorstr=$3
outpath=$4

IFS='+' read -r -a prior <<< "$priorstr"

sample-pop-params $popmodel -n $numsamps -p ${prior[@]} -o $outpath -v
