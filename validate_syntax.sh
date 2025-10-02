#!/bin/bash
# Quick syntax validation for Nextflow pipeline

echo "Validating Nextflow pipeline syntax..."

cd /data/daniel/genomedenovo_nf

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf_gp_env

# Try to parse the pipeline (this will fail early if there are syntax errors)
echo ""
echo "Running: nextflow run main.nf -config nextflow.config --help"
echo ""

nextflow run main.nf -config nextflow.config --help

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Pipeline syntax is valid!"
    echo ""
    echo "You can now run the actual pipeline with:"
    echo ""
    echo "  micromamba activate nf_gp_env"
    echo "  ./run.sh"
    echo ""
else
    echo ""
    echo "✗ Pipeline has syntax errors. Check the output above."
    exit 1
fi
