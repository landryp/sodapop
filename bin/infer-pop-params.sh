#!/bin/bash

popparams=$1
likestr=$2
colstr=$3
likesamps=$4
popmodelstr=$5
priorstr=$6
popsamps=$7
bhparams=$8
selectfunc=$9
selectpriorstr=${10}
selectsamps=${11}
numpost=${12}
numwalkers=${13}
numburnin=${14}
outpath=${15}
batch=${16}

IFS=',' read -r -a like <<< "$likestr"
IFS=',' read -r -a col <<< "$colstr"
IFS=',' read -r -a popmod <<< "$popmodelstr"
IFS='=' read -r -a prior <<< "$priorstr"
IFS=',' read -r -a bh <<< "$bhparams"
IFS='+' read -r -a selprior <<< "$selectpriorstr"

infer-pop-params $popparams ${like[@]} -c ${col[@]} -l $likesamps -p ${popmod[@]} -P ${prior[@]} -n $popsamps -B ${bh[@]} -f $selectfunc -S ${selprior[@]} -s $selectsamps -t $numpost -w $numwalkers -b $numburnin -o $outpath --batch $batch -v
