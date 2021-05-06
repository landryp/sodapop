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

infer-pop-params $popparams $likestr -c $colstr -l $likesamps -p $popmodelstr -P $priorstr -n $popsamps -B $bhparams -f $selectfunc -S $selectpriorstr -s $selectsamps -t $numpost -w $numwalkers -b $nburnin -o $outpath --batch $batch -v
