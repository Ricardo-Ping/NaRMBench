#!/bin/bash
module load apps/samtools/1.16.1
module load apps/bedtools/2.30.0

exptdir=$1
ref=$2
numcore=$3
# exptdir="/public/xump/202409Nanopore_tools/05_models/m5C/retrain/DInoPORE/01_feature/IVET/IVET_mod"
# ref="/public/xump/202409Nanopore_tools/05_models/m5C/retrain/DInoPORE/01_feature/IVET/IVET_reference.fa"
# numcore=24

codedir="/public/xump/202409Nanopore_tools/05_models/m5C/retrain/DInoPORE/code"
export SAM2TSV="java -jar ${codedir}/misc/sam2tsv.jar"
export PICARD="java -jar ${codedir}/misc/picard.jar"
export graphmap2="${codedir}/misc/graphmap2"

expt=$(basename ${exptdir})
fqdir=${exptdir}/out_fastq_bam
npdir=${exptdir}/out_nanopolish
dict=$(echo $ref | sed 's/.fasta//g' | sed 's/.fa//g').dict

echo =========================================================================================
echo "S2 Start: $(date)"

#Part 1 - sam2tsv
#Step 1 - Create Sequence Dictionary
if [ ! -f "$dict" ]; then
	echo Sequence dictionary not found. Creating sequence dictionary.
	$PICARD CreateSequenceDictionary R=$ref O=$dict
else
	echo Sequence dictionary found.
fi

#Step 2 - Run sam2tsv and process tsv file
sh ${codedir}/s2.Sam2tsv_processtsv.sh $fqdir $expt $exptdir $ref $numcore

#Part 2 - process raw nanopolish file
sh ${codedir}/s2.Combine_raw_nnpl.sh $npdir $expt $exptdir $numcore

echo -e S2 End: $(date) "\n"
echo =========================================================================================
