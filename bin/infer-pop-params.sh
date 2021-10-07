#!/bin/bash

popparams=$1
likestr=$2
classtr=$3
colstr=$4
likesamps=$5
popmodelstr=$6
priorstr=$7
distprior=$8
bhparams=$9
selectfunc=${10}
selectpriorstr=${11}
selectsamps=${12}
numpost=${13}
numwalkers=${14}
numburnin=${15}
outpath=${16}
batch=${17}
fixedstr=${18}

IFS=',' read -r -a like <<< "$likestr"
IFS=',' read -r -a clas <<< "$classtr"
IFS=',' read -r -a col <<< "$colstr"
IFS=',' read -r -a popmod <<< "$popmodelstr"
IFS='=' read -r -a prior <<< "$priorstr"
IFS=',' read -r -a bh <<< "$bhparams"
IFS='+' read -r -a selprior <<< "$selectpriorstr"
IFS='+' read -r -a fixed <<< "$fixedstr"

infer-pop-params $popparams ${like[@]} -C ${clas[@]} -c ${col[@]} -l $likesamps -p ${popmod[@]} -P ${prior[@]} -F ${fixed[@]} -D $distprior -B ${bh[@]} -f $selectfunc -S ${selprior[@]} -s $selectsamps -t $numpost -w $numwalkers -b $numburnin -o $outpath --batch $batch -v
