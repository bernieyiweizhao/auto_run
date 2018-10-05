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
  sleep 5
done

#set environment variables
export PATH=$PATH:/brcwork/bioinf/tools/bcl2fastq-2.20/bin
export PATH=$PATH:/brcwork/sequence/cellranger-2.2.0
export SGE_CLUSTER_NAME=brclogin1.cm.cluster

export PATH=$PATH:/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/bin #contains pandoc
Rscript="/brcwork/bioinf/tools/R/R-3.4.3/bin/Rscript"

#run mkfastq
col_names=$(grep -A 1 '\[Data]' $samplesheet|tail -1)
first_col=$(echo $col_names| awk -F ',' '{print $1}')

if [ $first_col = Lane ]; then
  echo "cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100 --lanes=1,2,3,4"
  echo "cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100 --lanes=1,2,3,4"|qsub -N "mkfastq_$mkfastqID"  -l h_vmem=6G,mem_free=6G
else
  echo "cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100"
  echo "cellranger mkfastq --id=$mkfastqID --run=$bcl --samplesheet=$samplesheet --jobmode=sge --maxjobs=100"|qsub -N "mkfastq_$mkfastqID"  -l h_vmem=6G,mem_free=6G
fi

#check for mkfastq completion
mkfastq_complete="$mkfastqID/_vdrkill"
while true; do
  if [ -f $mkfastq_complete ]; then
    echo "[$(date)] mkfastq complete"
    sleep 30
    break
  else
    echo "[$(date)] Running mkfastq"
  fi
  sleep 5
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
  countGenome=$(echo $line|awk -F ',' '{print $9}')
  curr=$(pwd)
  fastq_path="$curr/$mkfastqID/outs/fastq_path/"
  countID=$fastq
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
  output="$curr/$progress/$countID.output.txt"
  error="$curr/$progress/$countID.error.txt"
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"|qsub -N "count_$countID" -l h_vmem=6G,mem_free=6G -o $output -e $error
  
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
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
  output="$curr/$progress/$countID.output.txt"
  error="$curr/$progress/$countID.error.txt"
  echo "cellranger count --id=$countID --transcriptome=$countGenome --fastqs=$fastq_path --jobmode=sge --maxjobs=100"|qsub -N "count_$countID" -l h_vmem=6G,mem_free=6G -o $output -e $error
  count_outputs[$countID]="running"
done <<< "$lines"
fi

#run seurat
while [ ${#count_outputs[@]} -gt 0 ]; do
  for count_output in "${!count_outputs[@]}"; do
    echo "$count_output ${test[$key]}"
    if [ -f "$bcl/$count_output/_vdrkill" ]; then
      barcodes="$bcl/$count_output/outs/filtered_gene_bc_matrices/*/barcodes.tsv"
      num_cells=$(wc -l $barcodes|awk '{print $1}')
      mem=$(( $num_cells/1000 + 10 )) # 1gb for every thousand cells (floored) plus 10gb
      
      #submit standard seurat script
      command="$Rscript $seurat_dir/Seurat_main.R $bcl/$count_output $bcl/$count_output/outs $seurat_dir"
      output="$bcl/$progress/seurat_${count_output}.output.txt"
      error="$bcl/$progress/seurat_${count_output}.error.txt"
      echo $output
      echo $error 
      echo "$command"
      echo "$command"|qsub -N "seurat_${count_output}" -l h_vmem=${mem}G,mem_free=${mem}G -o $output -e $error
     
      # remove the current count output name from array
      unset count_outputs[$count_output]
    fi
  done
  sleep 10
done

echo "DONE"

