#!/bin/bash

#place this script in the output directory of sequencing run

mkfastqID=$1 #id for cellranger mkfastq output
samplesheet=$2 #name of the sample sheet file

#other global variable for paths
bcl=$(pwd)
RTAcomplete="RTAComplete.txt"
seurat_dir="/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/seurat_standard"

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
export PATH=/brcwork/bioinf/tools/bcl2fastq-2.20/bin:$PATH
export PATH=/brcwork/sequence/cellranger-2.2.0:$PATH
export SGE_CLUSTER_NAME=brclogin1.cm.cluster

export PATH=$PATH:/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/bin #contains pandoc
Rscript="/brcwork/bioinf/tools/R/R-3.4.3/bin/Rscript"

#check samplesheet format
col_names=$(grep -A 1 '\[Data]' $samplesheet|tail -1)
first_col=$(echo $col_names| awk -F ',' '{print $1}')

#run mkfastq
curr_mkfastq_job=-1
attempt=0

while true; do
  
  #start new attempt of mkfastq
  attempt=$(( $attempt+1 ))
  mkfastqID="${1}_$attempt"  
  
  if [ $first_col = Lane ]; then
    command="cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --localmem=64 --localcores=8 --lanes=1,2,3,4"
  else
    command="cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --localmem=64 --localcores=8"
  fi
  
  echo "[$(date)] $command"
  curr_mkfastq_job=$(echo "$command" |qsub -N "mkfastq_$mkfastqID"  -l h_vmem=8G,mem_free=64G -pe ncpus 8|awk '{print $3}')
  echo "[$(date)] Submitted mkfastq attempt $attempt; job id = $curr_mkfastq_job"

  #keep checking if current attempt is complete/stuck (mkfastq would sometimes get stuck in the BCL2FASTQ step; a re-run would usually work)
  while [ $(qstat|awk '{print $1}'|grep "^$curr_mkfastq_job"|wc -l) -gt 0 ]; do
    
    #kill job if stuck (if there are warning messeages in BCL2FASTQ step)
    err_file="$bcl/$mkfastqID/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-*/_stderr"
    if [ -f $err_file ]; then
      num_warning=$(grep "WARNING:" $err_file|wc -l)
      if [ $num_warning -gt 0 ]; then
        echo "[$(date)] mkfastq attempt $attempt failed"
        qdel $curr_mkfastq_job
        break        
      fi 
    fi
    
    echo "[$(date)] mkfastq attempt $attempt running"
    sleep 3m
  done
  
  #check if mkfastq is succesful
  if [ -f "$bcl/$mkfastqID/_vdrkill" ]; then
    echo "[$(date)] mkfastq attempt $attempt complete"
    sleep 30
    break
  fi

done

#directory of count job output
progress="count_progress"
mkdir $progress

#store count output names in an associated array (treat as a set)
declare -A count_outputs

#run count
if [ $first_col = Lane ]; then

  line=$(grep -A 2 '\[Data]' $samplesheet|tail -1|tr -d '\r')
  fastq=$(echo $line|awk -F ',' '{print $2}')
  project=$(echo $line|awk -F ',' '{print $8}')
  countGenome=$(echo $line|awk -F ',' '{print $9}')
  curr=$(pwd)
  fastq_path="$curr/$mkfastqID/outs/fastq_path/"
  countID="${project}__${fastq}"
  command="cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
  echo "[$(date)] $command"
  output="$curr/$progress/$countID.output.txt"
  error="$curr/$progress/$countID.error.txt"
  echo "$command"|qsub -N "count_$countID" -l h_vmem=6G,mem_free=6G -o $output -e $error
  
  count_outputs[$countID]="running"
else

  lines=$(grep -A $(wc -l $samplesheet|awk '{print $1}') '\[Data]' $samplesheet|tail -n+3)
  
  while read line; do
    line=$(echo "$line"|tr -d '\r')
    fastq=$(echo $line|awk -F ',' '{print $1}')
    project=$(echo $line|awk -F ',' '{print $7}')
    countGenome=$(echo $line|awk -F ',' '{print $8}')
    curr=$(pwd)
    fastq_path="$curr/$mkfastqID/outs/fastq_path/$project/$fastq"
    countID="${project}__${fastq}"
    command="cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
    echo "[$(date)] $command"
    output="$curr/$progress/$countID.output.txt"
    error="$curr/$progress/$countID.error.txt"
    echo "$command"|qsub -N "count_$countID" -l h_vmem=6G,mem_free=6G -o $output -e $error
    count_outputs[$countID]="running"
  done <<< "$lines"
fi

#run seurat
#
curr_seurat_job=-1

while [ ${#count_outputs[@]} -gt 0 ]; do
  
  #check if any count job have finished running
  echo "[$(date)] Count job running: $(echo ${!count_outputs[@]} |awk 'BEGIN {OFS=", "}{$1=$1}1')"
  for count_output in "${!count_outputs[@]}"; do
    if [ -f "$bcl/$count_output/_vdrkill" ]; then
      echo "[$(date)] Count job for $count_output complete"

      barcodes="$bcl/$count_output/outs/filtered_gene_bc_matrices/*/barcodes.tsv"
      num_cells=$(wc -l $barcodes|awk '{print $1}')
      #mem=$(( $num_cells/1000 + 10 )) # 1gb for every thousand cells (floored) plus 10gb
      mem=30 #seems to work for <10k cells
      
      #submit standard seurat script for the finished count output
      command="$Rscript $seurat_dir/Seurat_main.R $bcl/$count_output $bcl/$count_output/outs $seurat_dir"
      output="$bcl/$progress/seurat_${count_output}.output.txt"
      error="$bcl/$progress/seurat_${count_output}.error.txt"
      echo "[$(date)] $command"
      curr_seurat_job=$(echo "$command"|qsub -N "seurat_${count_output}" -l h_vmem=${mem}G,mem_free=${mem}G -o $output -e $error|awk '{print $3}')
      echo "[$(date)] Submitted seurat analysis job for $count_output; job id = $curr_seurat_job"

      # remove the finished count output name from array
      unset count_outputs[$count_output]
      
      #wait the current seurat job to finish running
      #sleep 5 #in case job does not show up immediately
      #while [ $(qstat|awk '{print $1}'|grep $curr_seurat_job|wc -l) -gt 0 ]; do 
      #  sleep 30
      #done
    fi
  done
  sleep 10m
done

echo "DONE"

