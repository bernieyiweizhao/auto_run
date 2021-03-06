library(rmarkdown)

Sys.setenv("RSTUDIO_PANDOC"="/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/bin/")

# Set file directories
args <- commandArgs(trailingOnly = T)
args[1] <- normalizePath(args[1])
args[2] <- normalizePath(args[2])

seurat_dir <- "/brcwork/sequence/10x_data/BernieWorkingDirectory/auto_run/seurat_standard"
input_dir <- args[1] #"/brcwork/sequence/10x_data/BernieWorkingDirectory/standard_script/pbmc3k"
output_dir <- args[2] #"/brcwork/sequence/10x_data/BernieWorkingDirectory/standard_script/pbmc3k"

output_dir <- sprintf("%s/seurat_analysis_results", output_dir)
dir.create(output_dir)
output_file <- sprintf("%s/seurat_analysis_results.html", output_dir)
rmd_path <- sprintf("%s/Seurat_main.Rmd", seurat_dir)

intermediates_dir <- output_dir #for .md files

rmarkdown::render(input = rmd_path, output_dir = output_dir, output_file = output_file, intermediates_dir = intermediates_dir, clean = T)


