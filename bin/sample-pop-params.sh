#!/bin/bash

popmodel=$1
numsamps=$2
priorstr=$3
fixedstr=$4
outpath=$5

IFS='=' read -r -a prior <<< "$priorstr"
IFS='+' read -r -a fixed <<< "$fixedstr"

sample-pop-params $popmodel -n $numsamps -p ${prior[@]} -F ${fixed[@]} -o $outpath -v
