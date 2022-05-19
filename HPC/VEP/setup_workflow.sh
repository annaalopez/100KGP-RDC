#!usr/bin/env bash

set -e
endDir=$1
mkdir -p ${endDir}/{logs}

workflowDir="/gel_data_resources/workflows/BRS_functional_annotation_workflow/v1.1/"

cp $workflowDir/submit_workflow.sh \
$workflowDir/cromwell.conf \
$workflowDir/workfflow_options.json \
$workflowDir/workflow_inputs_b37.json \
$workflowDir/workflow_inputs_b38.json \
$workflowDir/workflow_inputs_b38_commercial.json \
$workflowDir/BRS_functional_annotation_workflow.wdl \
${endDir}/