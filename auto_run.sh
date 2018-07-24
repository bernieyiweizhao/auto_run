#!/bin/bash

#place this script in the output directory of sequencing run

mkfastqID=$1 #id for cellranger mkfastq output
samplesheet=$2 #name of the sample sheet file

bcl=$(pwd)
RTAcomplete="RTAComplete.txt"

#check for sequencing run completion every 30 minutes
while true; do
  if [ -f $RTAcomplete ]; then
    echo "[$(date)] $RTAcomplete found"
    break
  else
    echo "[$(date)] Still Running"
  fi
  sleep 30m
done

#set environment variables
export PATH=$PATH:/brcwork/bioinf/tools/bcl2fastq-2.20/bin
export PATH=$PATH:/brcwork/sequence/cellranger-2.1.1
export SGE_CLUSTER_NAME=brclogin1.cm.cluster

#run mkfastq
first_col=$(grep -A 1 '\[Data]' $samplesheet|tail -1|awk -F ',' '{print $1}')

if [ $fist_col = Lane ]; then   
  cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100 --lanes=1,2,3,4
else
  cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100
fi

#directory of count job output
progress="count_progress"
mkdir $progress

#run count
grep -A $(wc -l $samplesheet|awk '{print $1}') '\[Data]' $samplesheet|tail -n+3|while read line; do
  line=$(echo "$line"|tr -d '\r')
  fastq=$(echo $line|awk -F ',' '{print $1}')
  project=$(echo $line|awk -F ',' '{print $7}')
  countGenome=$(echo $line|awk -F ',' '{print $8}')
  curr=$(pwd)
  fastq_path="$curr/$mkfastqID/outs/fastq_path/$project/$fastq"
  countID="${project}__${fastq}"
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
  output="$curr/$progress/$countID.output.txt"
  error="$curr/$progress/$countID.error.txt"
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"|qsub -N $countID -l h_vmem=6G,mem_free=6G -o $output -e $error
  if [ $first_col = Lane ]; then
    break
  fi
done
