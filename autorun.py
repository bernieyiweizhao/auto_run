import subprocess
import sys
import os
import pandas as pd
from io import StringIO
import numpy as np
import gzip
import time
import glob
import pathlib
import shutil
import re
##################

def print_flush(output):
  print(output, flush = True)
  pass

def get_time():
  T = time.localtime(time.time())
  return("[{}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}]".format(*T[:6]))

def set_cellranger_path(input):
  if input["workflow"] in ["rna", "vdj"]:
    if input["chemistry"] == "v2":
      cellranger="/brcwork/sequence/archived_cell_ranger_versions/cellranger-2.2.0/cellranger"
    else: #v3
      cellranger="/brcwork/sequence/cellranger_v3/cellranger-3.0.2/cellranger"
  else: #workflow = atac
    cellranger="/brcwork/sequence/cellranger_atac/cellranger-atac-1.1.0/cellranger-atac"
  return cellranger

def is_job_active(job_id):
  ps = subprocess.Popen(["qstat"], stdout = subprocess.PIPE)
  job_info = ps.communicate()[0].decode("utf-8")
  job_info =  pd.read_csv(StringIO(job_info), delim_whitespace=True)
  return(job_id in set(job_info.loc[:,"job-ID"]))

def extract_sample_info(samplesheet_path):
  samplesheet = pd.read_csv(samplesheet_path)
  data_line_num = samplesheet.index[samplesheet.iloc[:,0] == "[Data]"][0]
  sample_info = pd.read_csv(samplesheet_path,header = data_line_num + 2, dtype = "object")
  sample_info.fillna("", inplace=True)
  if "Lane" in list(sample_info):
    sample_info = sample_info.loc[[0]]
  for i in range(sample_info.shape[0]):
    if sample_info.loc[i,"CountID"] == "":
      sample_info.loc[i,"CountID"] = "{}__{}".format(sample_info.loc[i,"Sample_Project"], sample_info.loc[i,"Sample_Name"])
  return(sample_info)

def check_input():
  return

#get fastq path for each sample
def get_path(input, sample_info):
  output = ""
  for i in range(sample_info.shape[0]):
    project = sample_info.loc[i,"Sample_Project"]
    sample = sample_info.loc[i,"Sample_ID"]
    fastq_path = os.path.abspath("{}/{}".format(input["fastq_path"], project))
    output += "{}: {}\n".format(sample, fastq_path)
  return output

def run_md5sum(Path, md5_file):
  f = open(md5_file, 'a')
  for file in pathlib.Path(Path).glob("**/*"):
    file = str(file)
    if not os.path.isdir(file):
      ps = subprocess.Popen(["md5sum",file], stdout = subprocess.PIPE)
      f.write(ps.communicate()[0].decode("utf-8"))
  f.close()

def check_sequencing():

  print_flush("Start")
  while not os.path.isfile("RTAComplete.txt"):
    print_flush(get_time() + " Still Running")
    time.sleep(600)
  return

def is_mkfastq_complete(mkfastqID):
  return os.path.isfile("{}/_vdrkill".format(mkfastqID))

def is_mkfastq_error(mkfastqID):
  err_file = glob.glob("{}/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-*/_stderr".format(mkfastqID))
  if len(err_file) == 0:
    return False
  err = open(err_file[0], "r")
  for line in err:
    if ("WARNING:" in line) or ("Backtrace:" in line):
      err.close()
      return True
  err.close()
  return False

def run_mkfastq(input,sample_info):
  mkfastqID = "{}".format(input["mkfastqID"])
  if is_mkfastq_complete(mkfastqID):
    print_flush(get_time() + " mkfastq complete")
    return True

  attempt = 0
  max_attampt=5
  
  for attempt in range(1, max_attampt + 1):
    command = [input["cellranger"], "mkfastq", "--id={}".format(mkfastqID), "--run=.", "--samplesheet={}".format(input["samplesheet"]), "--localmem=36", "--localcores=6"]
    if "Lane" in list(sample_info):
      command = command + ["--lanes=1,2,3,4"] 
    print_flush(command)
    
    job_name = "mkfastq_" + mkfastqID
    output = "{}/{}.output.txt".format(input["log"], job_name)
    error = "{}/{}.error.txt".format(input["log"], job_name)
    ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
    ps = subprocess.Popen(["qsub","-N", job_name, "-l", "h_vmem=7G", "-pe", "ncpus", "6", "-o", output, "-e", error], stdin = ps.stdout, stdout = subprocess.PIPE)
    job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
    print_flush("job id = {}".format(job_id))
    while True:
      if is_mkfastq_complete(mkfastqID):
        print_flush(get_time() + " mkfastq complete")
        run_md5sum(mkfastqID + "/outs/fastq_path/", input["md5_file"])
        return True
      if is_mkfastq_error(mkfastqID):
        print_flush(get_time() + " mkfastq attempt {} failed".format(attempt))
        print_flush(job_id)
        subprocess.call(["qdel",job_id])
        os.rename(mkfastqID, "{}_{}".format(mkfastqID, attempt))
        break
      print_flush(get_time() + " mkfastq attempt {} running".format(attempt))
      time.sleep(60)

  print_flush(get_time() + " mkfastq reaches max number of attempts")
  return False

def is_trimming_complete(input):
  return os.path.isfile("{}/{}".format(input["trim_fastqID"], "_complete"))

def run_trimming(input, sample_info):
  cutadapt_jobs = set()
  complete = "_complete"
  if is_trimming_complete(input):
    return cutadapt_jobs

  fastq_path = "{}/{}/{}".format(input["mkfastqID"], "outs", "fastq_path")
  for i in range(sample_info.shape[0]):
    project = sample_info.loc[i,"Sample_Project"]
    sample = sample_info.loc[i,"Sample_ID"]
    trimmed_fastq_path = "{}/{}/{}".format(input["trim_fastqID"], project, sample)
    os.makedirs(trimmed_fastq_path, exist_ok = True)
    for r2_file in pathlib.Path(fastq_path).glob("**/{}*R2_001.fastq.gz".format(sample)):
      print_flush(r2_file)
      trimmed_r2_file = "{}/{}".format(trimmed_fastq_path, os.path.basename(r2_file))
      command = ["cutadapt", "--no-indels", "-e", "0.05", "-g", "^AAGCAGTGGTATCAACGCAGAGT", "-j", "4", "-o", trimmed_r2_file, r2_file]
      job_name = "cutadapt_" + os.path.basename(r2_file).split(".")[0]
      output = "{}/{}.output.txt".format(input["log"], job_name)
      error = "{}/{}.output.txt".format(input["log"], job_name)
      ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
      ps = subprocess.Popen(["qsub", "-N", job_name, "-l", "h_vmem=1G", "-pe", "ncpus", "4", "-o", output, "-e", error], stdin = ps.stdout, stdout = subprocess.PIPE)
      job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
      print_flush(job_id)
      cutadapt_jobs.add(job_id)
  for i in range(sample_info.shape[0]):
    project = sample_info.loc[i,"Sample_Project"]
    sample = sample_info.loc[i,"Sample_ID"]
    trimmed_fastq_path = "{}/{}/{}".format(input["trim_fastqID"], project, sample)
    for r1_file in pathlib.Path(fastq_path).glob("**/{}*R1_001.fastq.gz".format(sample)):
      shutil.copy(r1_file, trimmed_fastq_path)
  return cutadapt_jobs

def check_trimming(input, cutadapt_jobs):
  while len(cutadapt_jobs) > 0:
    print_flush(get_time() + " cutadapt running")
    finished_jobs = set()
    for job_id in cutadapt_jobs:
      if not is_job_active(job_id):
        finished_jobs.add(job_id)
    time.sleep(60)
    cutadapt_jobs -= finished_jobs

  if not is_trimming_complete(input):
    f = open("{}/{}".format(input["trim_fastqID"], "_complete"), 'w')
    f.write(get_time() + "\n")
    f.close()
  print_flush(get_time() + " cutadapt complete")
  run_md5sum(input["trim_fastqID"], input["md5_file"]) 
  return

def run_count(input, sample_info):
  workflow = input["workflow"]
  cellranger = input["cellranger"]
  count_jobs = {}
  for i in range(sample_info.shape[0]):
    if sample_info.loc[i,"Run_count"] == "0":
      continue
    project = sample_info.loc[i,"Sample_Project"]
    sample = sample_info.loc[i,"Sample_ID"]
    reference = sample_info.loc[i,"Reference"]
    library_file = sample_info.loc[i,"Library_file"]
    feature_file = sample_info.loc[i,"Feature_Barcode_file"]
    fastq_path = "{}/{}".format(input["fastq_path"], project)
    added_fastq_path = sample_info.loc[i,"Additional_fastq"]
    Count_ID = sample_info.loc[i,"CountID"]

    if len(library_file) > 0 and len(feature_file) > 0:
      # Cite-seq
      command = [cellranger, "count", "--id=" + Count_ID, "--transcriptome=" + reference, "--libraries=" + library_file,"--feature-ref=" + feature_file, "--jobmode=sge", "--maxjobs=100"]
    else:
      if len(added_fastq_path) > 0:
        fastq_path = [fastq_path] + added_fastq_path.replace(" ","").split(";")
        fastq_path = ",".join(fastq_path)
      command = {
            "rna": [cellranger, "count", "--id=" + Count_ID, "--transcriptome=" + reference, "--fastqs=" + fastq_path, "--sample=" + sample, "--jobmode=sge", "--maxjobs=100"],
            "vdj": [cellranger, "vdj", "--id=" + Count_ID, "--reference=" + reference, "--fastqs=" + fastq_path, "--sample=" + sample, "--jobmode=sge", "--maxjobs=100"],
            "atac": [cellranger, "count", "--id=" + Count_ID, "--reference=" + reference, "--fastqs=" + fastq_path, "--sample=" + sample, "--jobmode=sge", "--maxjobs=100"]
            }[workflow]
    
    job_name = "count_" + Count_ID
    output = "{}/{}.output.txt".format(input["log"], job_name)
    error = "{}/{}.error.txt".format(input["log"], job_name)
    ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
    ps = subprocess.Popen(["qsub", "-N", job_name, "-l", "h_vmem=8G", "-e", error, "-o", output], stdin = ps.stdout, stdout = subprocess.PIPE)
    job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
    print_flush(job_id)
    print_flush(command)
    count_jobs[sample] = job_id
  return count_jobs

def v3_to_v2(Count_ID):
  mat_path_v3 = Count_ID + "/outs/filtered_feature_bc_matrix/"
  mat_path = Count_ID + "/outs/filtered_gene_bc_matrices/genome/"
  os.makedirs(mat_path, exist_ok = True)
  barcodes = pd.read_csv(mat_path_v3 + "barcodes.tsv.gz", compression = "gzip", header = None, sep = "\t")
  barcodes.to_csv(mat_path + "barcodes.tsv", sep = "\t", header = None, index = False)
  genes = pd.read_csv(mat_path_v3 + "features.tsv.gz", compression = "gzip", header = None, sep = "\t")
  genes = genes.iloc[:,:2]
  genes.to_csv(mat_path + "genes.tsv", sep = "\t", header = None, index = False)
  f_mat_gz = gzip.open(mat_path_v3 + "matrix.mtx.gz")
  output = f_mat_gz.read()
  f_mat_gz.close()
  f_mat = open(mat_path + "matrix.mtx", 'w')
  f_mat.write(output.decode("utf-8"))
  f_mat.close()
  return

def rename_summary(Count_ID):
  os.rename("{}/outs/web_summary.html".format(Count_ID), "{}/outs/{}_web_summary.html".format(Count_ID, Count_ID))
  loupe = str(next(pathlib.Path(".").glob("{}/outs/*loupe".format(Count_ID))))
  os.rename(loupe, re.sub('outs/.loupe', "outs/{}".format(Count_ID), loupe))
  pass

def run_seurat(input, Count_ID):
  command = ["Rscript", input["seurat_dir"] + "/Seurat_main.R", Count_ID, Count_ID + "/outs"]
  job_name = "seurat_" + Count_ID
  output = "{}/{}.output.txt".format(input["log"], job_name)
  error = "{}/{}.error.txt".format(input["log"], job_name)
  ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
  ps = subprocess.Popen(["qsub","-N",job_name, "-l", "h_vmem=30G", "-o", output, "-e", error], stdin = ps.stdout, stdout = subprocess.PIPE)
  job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
  return job_id

def run_monocle(input, Count_ID):
  command = ["Rscript", input["monocle_dir"] + "/Monocle_main.R", Count_ID, Count_ID + "/outs"]
  job_name = "monocle_" + Count_ID
  output = "{}/{}.output.txt".format(input["log"], job_name)
  error = "{}/{}.error.txt".format(input["log"], job_name)
  ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
  ps = subprocess.Popen(["qsub","-N",job_name, "-l", "h_vmem=12G", "-pe", "ncpus", "4", "-o", output, "-e", error], stdin = ps.stdout, stdout = subprocess.PIPE)
  job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
  return job_id

def run_HTODemux(input, Count_ID, hashtag_path):
  command = ["Rscript", input["HTOdemux_dir"] + "/HTO_demux.R", Count_ID, Count_ID + "/outs", hashtag_path, Count_ID]
  job_name = "HTOdemux_" + Count_ID
  output = "{}/{}.output.txt".format(input["log"], job_name)
  error = "{}/{}.error.txt".format(input["log"], job_name)
  ps = subprocess.Popen(["echo"] + command, stdout = subprocess.PIPE)
  ps = subprocess.Popen(["qsub","-N",job_name, "-l", "h_vmem=8G", "-o", output, "-e", error], stdin = ps.stdout, stdout = subprocess.PIPE)
  job_id = ps.communicate()[0].decode("utf-8").strip().split(" ")[2]
  return job_id

def run_secondary(input, sample_info, count_jobs):
  seurat_jobs = {}
  monocle_jobs = {}
  HTODemux_jobs = {}
  while True:
    progress = ""
    for i in range(sample_info.shape[0]):
      if sample_info.loc[i,"Run_count"] == "0":
        continue

      project = sample_info.loc[i,"Sample_Project"]
      sample = sample_info.loc[i,"Sample_ID"]
      Count_ID = sample_info.loc[i,"CountID"]
      if (sample in count_jobs):
        if is_job_active(count_jobs[sample]):
          progress += "\t{}: running count\n".format(sample)
        else:
          del count_jobs[sample]
          if input["workflow"] == "rna":
            if len(sample_info.loc[i,"Library_file"]) > 0: #Cite-seq
              v3_to_v2(Count_ID)
              HTODemux_jobs[sample] = run_HTODemux(input, Count_ID, sample_info.loc[i,"Feature_Barcode_file"])
            else:
              if input["chemistry"] == "v3":
                v3_to_v2(Count_ID)
              seurat_jobs[sample] = run_seurat(input, Count_ID)
              monocle_jobs[sample] = run_monocle(input, Count_ID)
          
          rename_summary(Count_ID)
          run_md5sum(Count_ID + "/outs/", input["md5_file"])

      if (sample in HTODemux_jobs):
        if is_job_active(HTODemux_jobs[sample]):
          progress += "\t{}: running HTOdemux\n".format(sample)
        else:
          del HTODemux_jobs[sample]
          run_md5sum(Count_ID + "/outs/HTOdemux_results/", input["md5_file"])

      if (sample in seurat_jobs):
        if is_job_active(seurat_jobs[sample]):
          progress += "\t{}: running seurat\n".format(sample)
        else:
          del seurat_jobs[sample]
          run_md5sum(Count_ID + "/outs/seurat_analysis_results/", input["md5_file"])

      if (sample in monocle_jobs):
        if is_job_active(monocle_jobs[sample]):
          progress += "\t{}: running monocle\n".format(sample)
        else:
          del monocle_jobs[sample]
          run_md5sum(Count_ID + "/outs/monocle_analysis_results/", input["md5_file"])

    if len(count_jobs) + len(seurat_jobs) + len(monocle_jobs) + len(HTODemux_jobs) == 0:
      break
    print_flush(get_time() + "\n" + progress.strip("\n"))
    time.sleep(300)

  run_md5sum(CountID + "/outs/", input["md5_file"])
  print_flush(get_time() + " DONE")
  return sample_info

def main():
  input={}
  action = sys.argv[1]
  input["mkfastqID"]=sys.argv[2]
  input["samplesheet"]=sys.argv[3]
  input["workflow"]=sys.argv[4]
  input["chemistry"]=sys.argv[5]
  print_flush(input)
  
  ### envrioment variables
  os.environ["PATH"] = "/brcwork/bioinf/tools/pigz-2.4/:/brcwork/bioinf/tools/R/R-3.4.3/bin/:/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/bin:/brcwork/bioinf/tools/bcl2fastq-2.20/bin:" + os.environ["PATH"]
  os.environ["SGE_CLUSTER_NAME"] = "brclogin1.cm.cluster"

  ### global variables

  input["cellranger"] = set_cellranger_path(input)
  input["seurat_dir"] = "/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/seurat_standard"
  input["monocle_dir"] = "/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/monocle_standard"
  input["HTOdemux_dir"] = "/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/HTODemux_standard/"
  input["log"]="auto_run_log"
  input["md5_file"] = "checksum.md5"
  input["trim_fastqID"] = "{}.trimmed".format(input["mkfastqID"])
  input["run_trim"] = input["workflow"] == "rna"
  if input["run_trim"]:
    input["fastq_path"] = input["trim_fastqID"]
  else:
    input["fastq_path"] = input["mkfastqID"] + "/outs/fastq_path"

  ### run

  sample_info = extract_sample_info(input["samplesheet"]) 
  print_flush(sample_info)
  if action == "run":
    print_flush(input)
    check_input()
    check_sequencing()
    os.makedirs(input["log"], exist_ok = True)  
    run_mkfastq(input, sample_info)
    if input["run_trim"]:
      cutadapt_jobs = run_trimming(input, sample_info)
      check_trimming(input,cutadapt_jobs) 
    count_jobs = run_count(input, sample_info)
    run_secondary(input, sample_info, count_jobs)
  elif action == "get_path":
    print_flush(get_path(input, sample_info))
  else:
    print_flush("Invalid action")
  return

main()

