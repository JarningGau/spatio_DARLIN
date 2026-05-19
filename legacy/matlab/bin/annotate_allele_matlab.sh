#!/bin/bash

CARLIN_dir=$1
input_dir=$2
output_dir=$3
sample=$4
cfg_type=$5
template=$6
read_cutoff_UMI_override=$7
read_cutoff_CB_override=$8

# # This is used to test the pipeline
# CARLIN_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN'
# input_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN/tests/Tigre_Carlin_test'
# output_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN/tests/Tigre_Carlin_test/output2'
# sample='test_Tigre_Carlin'
# cfg_type='BulkDNA_Tigre'
# template='Tigre'
# read_cutoff_CB_override=3
# read_cutoff_UMI_override=8

echo "Running interactive mode for CARLIN analysis"
command_str="my_CARLIN_pipeline('$sample','$cfg_type','$input_dir','$output_dir','$template','read_cutoff_UMI_override',$read_cutoff_UMI_override,'read_cutoff_CB_override',$read_cutoff_CB_override,'CARLIN_dir','$CARLIN_dir')"
echo $command_str
cur_dir=$(pwd)
cd $CARLIN_dir
matlab -nodisplay -nosplash -nodesktop -r "$command_str; exit"
cd $cur_dir
