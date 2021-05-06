#!/bin/bash

popmodel=$1
numsamps=$2
priorstr=$3
outpath=$4

sample-pop-params $popmodel -n $numsamps -p $priorstr -o $outpath -v
