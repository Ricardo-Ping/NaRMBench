#!/bin/bash

exptdir=$1
agggrp=$2
numcore=$3

codedir=$(cd `dirname $0` && pwd)

input=$agggrp.morefts.input_CNN_regression_modgen.RData

input2=$agggrp.morefts.input_CNN_regression_modgen_noclassref.RData

rpath=${codedir}/s6.Predict_test_data.R

rpath2=${codedir}/s6.Predict_test_data_noclassref.R

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S6 Start: $(date)"


if [ -f "$agggrp.morefts.input_CNN_regression_modgen_noclassref.RData" ]; then
    Rscript $rpath2 -i $input2 -t $numcore
else
    Rscript $rpath -i $input -t $numcore
fi
#Rscript $rpath -i $input -t $numcore
echo -e S6 End: $(date) "\n"
echo =========================================================================================