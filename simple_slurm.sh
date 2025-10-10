#!/bin/bash
#SBATCH --job-name=genome_assembly_10
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_output_%j.log
#SBATCH --error=slurm_error_%j.log

# Activate environment
source /home/benjamin/micromamba/etc/profile.d/micromamba.sh
micromamba activate nf_gp_env

# Define your 10 isolates with their file patterns
declare -a ISOLATES=(
    "SCA3582A1_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "SCA3582A9_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "O68_*_L2_{1,2}.fq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Novogene_2020/X204SC20061145-Z01-F006/raw_data/O68/"
    "SCA3582A17_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "SCA3582A25_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "SCA3582A33_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "SCA3582A41_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "SCA3582A49_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
    "R1_*_L2_{1,2}.fq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Novogene_2020/X204SC20061145-Z01-F006/raw_data/R1/"
    "SCA3582A57_R{1,2}_001.fastq.gz|/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Syngenta_2020/200911_A00306_0043_AHHCK7DSXY/"
)

# Simple for loop through all isolates
for i in "${!ISOLATES[@]}"; do
    echo "Processing isolate $((i+1))/10: ${ISOLATES[$i]}"
    
    # Split the pattern and path
    IFS='|' read -r PATTERN PATH <<< "${ISOLATES[$i]}"
    FULL_PATH="${PATH}${PATTERN}"
    
    echo "Running pipeline on: $FULL_PATH"
    
    # Run nextflow for this isolate
    nextflow run main.nf \
        -config nextflow.config \
        -profile slurm \
        --outdir "./results_isolate_$((i+1))" \
        --reads "$FULL_PATH" \
        --blast_query /legserv/Temp/Benjamin/Avrstb6_sequence_IPO323_revcomp.fa \
        -resume
    
    echo "Completed isolate $((i+1))"
    echo "---"
done

echo "All 10 isolates completed!"