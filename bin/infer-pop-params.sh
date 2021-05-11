#!/bin/bash

popparams=$1
likestr=$2
colstr=$3
likesamps=$4
popmodelstr=$5
priorstr=$6
bhparams=$7
selectfunc=$8
selectpriorstr=$9
selectsamps=${10}
numpost=${11}
numwalkers=${12}
numburnin=${13}
outpath=${14}
batch=${15}

IFS=',' read -r -a like <<< "$likestr"
IFS=',' read -r -a col <<< "$colstr"
IFS=',' read -r -a popmod <<< "$popmodelstr"
IFS='=' read -r -a prior <<< "$priorstr"
IFS=',' read -r -a bh <<< "$bhparams"
IFS='+' read -r -a selprior <<< "$selectpriorstr"

infer-pop-params $popparams ${like[@]} -c ${col[@]} -l $likesamps -p ${popmod[@]} -P ${prior[@]} -B ${bh[@]} -f $selectfunc -S ${selprior[@]} -s $selectsamps -t $numpost -w $numwalkers -b $numburnin -o $outpath --batch $batch -v
