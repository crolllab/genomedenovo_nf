#!/bin/bash
#SBATCH --job-name=test_one
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=test_one_%j.log

# Activate environment
source /home/benjamin/micromamba/etc/profile.d/micromamba.sh
micromamba activate nf_gp_env

echo "Starting single isolate test"

# Just test the first isolate
nextflow run main.nf \
    -config nextflow.config \
    -profile slurm \
    --outdir ./slurm_test \
    --reads "/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/SCA3582A1_R{1,2}_001.fastq.gz" \
    --blast_query /legserv/Temp/Benjamin/Avrstb6_sequence_IPO323_revcomp.fa

echo "Finished!"