#!/usr/bin/Rscript
# This script can be used to analyze a batch of sequences without having to use
# GUI or R directly.
# It checks to see if the required tools are installed (devtools, mlplasmids)

if(!"devtools" %in% rownames(installed.packages())) {
    print("Installing devtools...")
    install.packages("devtools", repos='http://cran.us.r-project.org')
}

if(!"Biostrings" %in% rownames(installed.packages())) {
    print("Installing Biostrings...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install(version = "3.14")
	BiocManager::install("biocLite")
	BiocManager::install()
	BiocManager::install("Biostrings")   
}
if(!"mlplasmids" %in% rownames(installed.packages())) {
  print("Installing mlplasmids; please be patient, as this involves downloading a large dataset...")

  devtools::install_git("https://gitlab.com/mmb-umcu/mlplasmids.git")
}
suppressMessages(library(mlplasmids))

usage =  "USAGE: Rscript run_mlplasmids.R ./path/to/assembly.fasta ./path/to/output.tab [prob_threshold] [species]"
# ENABLE command line arguments
args <- commandArgs(TRUE)


# check the required arguments
input_path <- args[1]
output_path <- args[2]
if(any(is.na(c(input_path, output_path)))){
    print(usage)
    stop()
}

# set the defaults
thresh <- ifelse(!is.na(args[3]), as.numeric(args[3]), .8)
species <- ifelse(!is.na(args[4]), args[4], "Escherichia coli")
full_output <- ifelse(!is.na(args[5]), args[5], FALSE)
print(paste("Threshold:", thresh))
print(paste("Species:", species))
print(paste("Full output:", full_output))

example_prediction <- plasmid_classification(path_input_file = input_path,  prob_threshold=thresh, species = species, full_output = full_output)
if (is.null(example_prediction)){
    stop("Issue with mlplasmids; please try running interactively")
}

write.table(x=example_prediction, file=output_path, row.names=F, sep="\t")
