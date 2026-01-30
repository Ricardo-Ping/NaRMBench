#!/bin/bash

exptdir=$1
numcore=$2
agggrp=$3
classref=$4

codedir=$(cd `dirname $0` && pwd)

input=$agggrp.Agg.morefts.10bin.inML.txt
output=$agggrp.morefts.input_CNN_regression_modgen.RData
output2=$agggrp.morefts.input_CNN_regression_modgen_noclassref.RData

rpath=${codedir}/s5.Preprocess_data_matrix_inputCNN.R
rpath2=${codedir}/s5.Preprocess_data_matrix_inputCNN_noclassref.R

cnndir=$(dirname $exptdir)/matrix_CNN

cd $cnndir

echo =========================================================================================
echo "S5 Start: $(date)"

# Check if classref is provided (not empty)
if [ -z "$classref" ]; then
  # If classref is not provided, run the second script
  Rscript $rpath2 -t $numcore -i $input -o $output2
else
  # If classref is provided, run the first script with classref
  Rscript $rpath -t $numcore -i $input -o $output -c $classref
fi

echo -e S5 End: $(date) "\n"
echo =========================================================================================
