 ################################################################################
import os
import yaml
import glob
import pandas as pd

################################################################################
#### Set up ####

# This is the master snakefile that will call all the other snakefiles

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "./config_AMF.yaml"


metadata = pd.read_csv(config['metadata_file_path'], usecols=['SampleID'])
samples_all = metadata['SampleID'].tolist()

sample_dir = "/home/amhorst/2021_tomato_rhizo/RNA/reads/hostrm/"

samples = [
    s for s in samples_all
    if os.path.exists(f"{sample_dir}{s}_R1_rmhost.fq.gz") and os.path.exists(f"{sample_dir}{s}_R2_rmhost.fq.gz")
]

print(f"Samples to be processed: {samples}")


################################################################################
#### Rule all ####

rule all:
    input:
        #check_sortmerna = expand("Checks/SortmeRNA/{sample}_sortmeRNA.done", sample=samples),
        check_repair = expand("Checks/SortmeRNA/{sample}_repair.done", sample=samples),
        check_salmon_index = "Checks/salmon_index/salmon_index.done",
        check_salmon_quant = expand("Checks/salmon_index/salmon_quant_{sample}.done", sample=samples),
        check_LongOrfs = "Checks/transdecoder/check_LongOrfs.done",
        check_Predict = "Checks/transdecoder/check_Predict.done",
        check_salmon_merge_tpm = "Checks/salmon_merge/salmon_merge_tpm.done",
        check_salmon_merge_counts = "Checks/salmon_merge/salmon_merge_counts.done",
        check_eggnog_mapper = "Checks/eggnog_mapper/check_eggnog_mapper.done",

################################################################################
#### Pipeline ####

include: "10_rnaseq_annotation.smk"

