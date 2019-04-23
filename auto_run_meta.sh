#!/bin/bash

#place this script in the output directory of sequencing run

### input args
mkfastqID=$1 #id for cellranger mkfastq output
samplesheet=$2 #name of the sample sheet file
workflow=$3 # rna/vdj/atac
chemistry=$4 # v2/v3

if [[ $workflow =~ (rna|vdj) ]]; then
  if [ $chemistry = v2 ]; then
    cellranger="/brcwork/sequence/archived_cell_ranger_versions/cellranger-2.2.0/cellranger"
    override_path=""
  else #v3
    cellranger="/brcwork/sequence/cellranger_v3/cellranger-3.0.2/cellranger"
    if [ $workflow = rna ]; then
      override_path="/brcwork/sequence/cellranger_v3/cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/sge_override.json"
    else #vdj
      override_path=""
    fi
  fi
else #workflow = atac
  #cellranger="/brcwork/sequence/cellranger_atac/cellranger-atac-1.0.1/cellranger-atac"
  #override_path="/brcwork/sequence/cellranger_atac/cellranger-atac-1.0.1/martian-cs/v3.1.0/jobmanagers/sge_override.json"
  cellranger="/brcwork/sequence/cellranger_atac/cellranger-atac-1.1.0/cellranger-atac"
  override_path=""

fi

### other global variables
RTAcomplete="RTAComplete.txt"
seurat_dir="/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/seurat_standard"
monocle_dir="/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/monocle_standard"
log="auto_run_log"

### check for sequencing run completion every 30 minutes
while true; do
  if [ -f $RTAcomplete ]; then
    echo "[$(date)] $RTAcomplete found"
    break
  else
    echo "[$(date)] Still Running"
  fi
  sleep 30m
done

### set environment variables
export PATH=/brcwork/bioinf/tools/bcl2fastq-2.20/bin:$PATH
export SGE_CLUSTER_NAME=brclogin1.cm.cluster
export PATH=/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/bin:$PATH #contains pandoc
export PATH=/brcwork/bioinf/tools/R/R-3.4.3/bin/:$PATH

### process samplesheet
dos2unix $samplesheet
col_names=($(grep -A 1 '\[Data]' $samplesheet| tail -1| awk 'BEGIN{FS=","}{$1=$1}1'))

data_start=$(grep -n '\[Data]' $samplesheet| cut -d ":" -f 1)
data_end=$(wc -l $samplesheet| cut -d " " -f 1)
data_nrow=$(($data_end - $data_start - 1))

if [ ${col_names[0]} = Lane ]; then
  data_rows=1
else
  data_rows=$(seq $data_nrow)
fi

declare -A data
for i in ${!col_names[@]}; do
  for j in $data_rows; do
    data[${col_names[$i]},$j]=$(tail -$data_nrow $samplesheet|cut -d "," -f $(($i+1)) | head -$j | tail -1)
    echo ${col_names[$i]},$j: ${data[${col_names[$i]},$j]}
  done
done

### log directory
mkdir "$log"

### run mkfastq
curr_mkfastq_job=-1
attempt=0
max_attempt=5

while true; do
  
  #start new attempt of mkfastq
  attempt=$(( $attempt+1 ))

  if [ $attempt -gt $max_attempt ]; then
    echo "[$(date)] Too many attempts"
    exit
  fi

  mkfastqID="${1}_$attempt"  
  if [ ${col_names[0]} = Lane ]; then
    lanes="--lanes=1,2,3,4"
  fi
  
  command="$cellranger mkfastq --id=$mkfastqID --run=. --samplesheet=$samplesheet --localmem=36 --localcores=6 $lanes"
  
  echo "[$(date)] $command"
  curr_mkfastq_job=$(echo "$command" |qsub -N "mkfastq_$mkfastqID"  -l h_vmem=8G,mem_free=36G -pe ncpus 6 -o "$log/mkfastq_${mkfastqID}.output.txt" -e "$log/mkfastq_${mkfastqID}.error.txt"  |awk '{print $3}')
  echo "[$(date)] Submitted mkfastq attempt $attempt; job id = $curr_mkfastq_job"

  #keep checking if current attempt is complete/stuck (mkfastq would sometimes get stuck in the BCL2FASTQ step; a re-run would usually work)
  while [ $(qstat|awk '{print $1}'|grep "^$curr_mkfastq_job"|wc -l) -gt 0 ]; do
    
    #kill job if stuck (if there are warning messeages in BCL2FASTQ step)
    err_file="$mkfastqID/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-*/_stderr"
    if [ -f $err_file ]; then
      num_warning=$(grep "WARNING:" $err_file|wc -l)
      if [ $num_warning -gt 0 ]; then
        echo "[$(date)] mkfastq attempt $attempt failed"
        qdel $curr_mkfastq_job
        break        
      fi 
    fi
    
    echo "[$(date)] mkfastq attempt $attempt running"
    sleep 1m
  done
  
  #check if mkfastq is succesful
  if [ -f "$mkfastqID/_vdrkill" ]; then
    echo "[$(date)] mkfastq attempt $attempt complete"
    sleep 10
    break
  fi

done

### md5sum for fastqs
find $mkfastqID/outs/fastq_path/ -type f -exec md5sum {} \; > checksum.md5

### store count progress in an associated array (treat as a set)
unset count_jobs; declare -A count_jobs
unset seurat_jobs; declare -A seurat_jobs
unset monocle_jobs; declare -A monocle_jobs
unset progress; declare -A progress

### run count
for i in $data_rows; do
  sample_id=${data[Sample_ID,$i]}
  project=${data[Sample_Project,$i]}
  reference=${data[Description,$i]}
  if [ ${col_names[0]} = Lane ]; then
    fastq_path="$mkfastqID/outs/fastq_path/$project"
  else
    fastq_path="$mkfastqID/outs/fastq_path/$project/$sample_id"
  fi

  countID="${project}__${sample_id}"
  
  if [ $workflow = "rna" ]; then
    command="$cellranger count --id=$countID --transcriptome=$reference --fastqs=$fastq_path --jobmode=sge --maxjobs=100 --override=$override_path"
  elif [ $workflow = "vdj" ]; then
    command="$cellranger vdj --id=$countID --reference=$reference --fastqs=$fastq_path --jobmode=sge --maxjobs=100"
  else #workflow = atac
    command="$cellranger count --id=$countID --reference=$reference --fastqs=$fastq_path --jobmode=sge --maxjobs=100 --override=$override_path"
  fi  
  
  echo "[$(date)] $command"
  output="$log/$countID.output.txt"
  error="$log/$countID.error.txt"
  job_id=$(echo "$command"|qsub -N "count_$countID" -l h_vmem=6G,mem_free=6G -o $output -e $error|awk '{print $3}')
  count_jobs[$countID]=$job_id
  progress[$countID]="running_count"
done


### check analysis step of each sample
while [ ${#progress[@]} -gt 0 ]; do
  
  # check if any job have finished running
  for sample in ${!progress[@]}; do
    
    if [ ${progress[$sample]} = "running_count" ]; then
      # check for count completion
      if [ -f "$sample/_vdrkill" ]; then
        unset count_jobs[$sample]
        echo "[$(date)] Count job for $sample complete"

        if [ $workflow = rna ]; then

          if [ $chemistry = v3 ]; then
            # restore to v2 output format for Seurat and Monocle
            mat_path_v3="$sample/outs/filtered_feature_bc_matrix"
            mat_path="$sample/outs/filtered_gene_bc_matrices"
            mkdir $mat_path
            mkdir $mat_path/genome
            zcat $mat_path_v3/barcodes.tsv.gz > $mat_path/genome/barcodes.tsv
            zcat $mat_path_v3/features.tsv.gz > $mat_path/genome/features.tsv
            awk 'BEGIN{OFS="\t"}{print $1,$2}' $mat_path/genome/features.tsv > $mat_path/genome/genes.tsv
            zcat $mat_path_v3/matrix.mtx.gz > $mat_path/genome/matrix.mtx 
          fi

          progress[$sample]="running_secondary"
          
          barcodes="$sample/outs/filtered_gene_bc_matrices/*/barcodes.tsv"
          num_cells=$(wc -l $barcodes|awk '{print $1}')
          mem_seurat=30 #seems to work for <10k cells
          mem_monocle=12 #need to change
          num_core_monocle=4 #need to change
          
          # run Seurat
          command="Rscript $seurat_dir/Seurat_main.R $sample $sample/outs"
          output="$log/seurat_${sample}.output.txt"
          error="$log/seurat_${sample}.error.txt"
          echo "[$(date)] $command"
          seurat_job_id=$(echo "$command"|qsub -N "seurat_${sample}" -l h_vmem=${mem_seurat}G,mem_free=${mem_seurat}G -o $output -e $error|awk '{print $3}')
          echo "[$(date)] Submitted seurat analysis job for $sample; job id = $seurat_job_id"
          seurat_jobs[$sample]=$seurat_job_id
          
          # run monocle
          command="Rscript $monocle_dir/Monocle_main.R $sample $sample/outs"
          output="$log/monocle_${sample}.output.txt"
          error="$log/monocle_${sample}.error.txt"
          echo "[$(date)] $command"
          monocle_job_id=$(echo "$command"|qsub -N "monocle_${sample}" -l h_vmem=${mem_monocle}G,mem_free=$(( $mem_monocle*$num_core_monocle ))G -pe ncpus $num_core_monocle -o $output -e $error|awk '{print $3}')
          echo "[$(date)] Submitted monocle analysis job for $sample; job id = $monocle_job_id"
          monocle_jobs[$sample]=$monocle_job_id 
        fi
      fi

    else # running secondary
      # check for seurat completion
      if [ -f "$sample/outs/seurat_analysis_results/seurat_analysis_results.html" ]; then
        if ! [ -z ${seurat_jobs[$sample]} ]; then
          unset seurat_jobs[$sample] 
          echo "[$(date)] Seurat job for $sample complete" 
        fi
      fi
      
      # check for monocle completion
      if [ -f "$sample/outs/monocle_analysis_results/monocle_analysis_results.html" ]; then
        if ! [ -z ${monocle_jobs[$sample]} ]; then
          unset monocle_jobs[$sample]
          echo "[$(date)] Monocle job for $sample complete"
        fi
      fi
    fi
    
    # check if all jobs have finished for this sample
    if [[ -z ${count_jobs[$sample]} && -z ${seurat_jobs[$sample]} && -z ${monocle_jobs[$sample]} ]]; then
      unset progress[$sample]
      
      # md5sum for outs
      find $sample/outs/ -type f -exec md5sum {} \; >> checksum.md5 

    fi
  done
  
  echo "[$(date)]"
  for sample in ${!progress[@]}; do
    job_ids=( ${count_jobs[$sample]} ${seurat_jobs[$sample]} ${monocle_jobs[$sample]} )
    job_ids=$(IFS=, ; echo "${job_ids[*]}")
    echo "    $sample: ${progress[$sample]}, job ids = $job_ids"
  done

  sleep 10m
  
done

echo "[$(date)] DONE"

