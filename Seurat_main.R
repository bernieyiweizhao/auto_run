library(rmarkdown)

# Set file directories
args <- commandArgs(trailingOnly = T)
args[1] <- normalizePath(args[1])

seurat_dir <- "C:/Users/Yiwei Zhao/Desktop/auto_run/" #can be changed
input_dir <- args[1] #"/brcwork/sequence/10x_data/BernieWorkingDirectory/standard_script/pbmc3k"
output_dir <- args[1] #"/brcwork/sequence/10x_data/BernieWorkingDirectory/standard_script/pbmc3k"

output_dir <- sprintf("%s/seurat_analysis_results", output_dir)
dir.create(output_dir)
output_file <- sprintf("%s/seurat_analysis_results.html", output_dir)
rmd_path <- sprintf("%s/Seurat_main.Rmd", seurat_dir)

rmarkdown::render(input = rmd_path, output_dir = output_dir, output_file = output_file)

