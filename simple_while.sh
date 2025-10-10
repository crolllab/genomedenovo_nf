#!/bin/bash
#SBATCH --job-name=genome_5_isolates
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_5isolates_%j.log
#SBATCH --error=slurm_5isolates_%j.err

# Activating the environment
source /home/benjamin/micromamba/etc/profile.d/micromamba.sh
micromamba activate nf_gp_env

# Listing first 5 isolates      
Isolate_1="SCA3582A1_R{1,2}_001.fastq.gz"
Isolate_2="SCA3582A9_R{1,2}_001.fastq.gz"
Isolate_3="O68_*_L2_{1,2}.fq.gz"
Isolate_4="SCA3582A17_R{1,2}_001.fastq.gz"
Isolate_5="SCA3582A25_R{1,2}_001.fastq.gz"

# The paths on legserv
Path_Syngenta="/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
Path_O68="/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Novogene_2020/X204SC20061145-Z01-F006/raw_data/O68/"
Path_R1="/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Novogene_2020/X204SC20061145-Z01-F006/raw_data/R1/"

# A simple while loop for first 5 isolates
COUNTER=1
while [ $COUNTER -le 5 ]; do
    echo "Processing isolate $COUNTER/5"
    if [ $COUNTER -eq 1 ]; then
        READS="${Path_Syngenta}${Isolate_1}"
    elif [ $COUNTER -eq 2 ]; then
        READS="${Path_Syngenta}${Isolate_2}"
    elif [ $COUNTER -eq 3 ]; then
        READS="${Path_O68}${Isolate_3}"
    elif [ $COUNTER -eq 4 ]; then
        READS="${Path_Syngenta}${Isolate_4}"
    elif [ $COUNTER -eq 5 ]; then
        READS="${Path_Syngenta}${Isolate_5}"
    fi
    
    echo "Running: $READS"
    
    nextflow run main.nf \
        -config nextflow.config \
        -profile slurm \
        --outdir "./results_$COUNTER" \
        --reads "$READS" \
        --blast_query /legserv/Temp/Benjamin/Avrstb6_sequence_IPO323_revcomp.fa
    
    echo "Completed isolate $COUNTER"
    echo "---"
    COUNTER=$((COUNTER + 1))
done

echo "Done!"