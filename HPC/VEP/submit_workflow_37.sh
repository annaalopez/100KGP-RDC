#BSUB -q inter
#BSUB -P bio
#BSUB -o logs/%J.stdout
#BSUB -e logs/%J.stderr
BSUB -cwd /re_gecip/enhanced_interpretation/AnnaLopez/annotation/
#BSUB -n 2
#BSUB -R "rusage[mem=36000] span[hosts=1"

module load bio/cromwell/51
mkdir -p logs
tstamp=tstamp=$(date +%Y_%m_%d_%Hh%M)

export TMPDIR = /re_gecip/enhanced_interpretation/AnnaLopez/annotation/

java -Dconfig.file=cromwell.conf -jar $CROMWELL_JAR \
run functional_annotation_workflow.wdl \
-i workflow_inputs_b37.json \
-o workflow_options.json \
-m metadata_"$tstamp".json