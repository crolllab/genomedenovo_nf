#!/bin/bash
# Test script for multi-sample scaffold input

set -euo pipefail

echo "========================================="
echo "Testing Multi-Sample Scaffold Processing"
echo "========================================="

# Parameters
ASSEMBLIES_DIR="/data/daniel/wagABC_blast/assemblies"
QUERY="./wagA_gene.fa"
OUTDIR="./nf_output_multisample"
WORK_DIR="/scratch/daniel_tmp_multisample"

# Clean up previous test
echo "Cleaning up previous test outputs..."
rm -rf "$OUTDIR" "$WORK_DIR"

# Activate environment
echo "Activating micromamba environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf_gp_env

# Count expected samples
echo ""
echo "Counting samples in $ASSEMBLIES_DIR..."
sample_count=$(find "$ASSEMBLIES_DIR" -maxdepth 2 -name "scaffolds.fasta" | wc -l)
echo "Found $sample_count samples with scaffolds.fasta files"

# Show first 10 samples
echo ""
echo "First 10 samples:"
find "$ASSEMBLIES_DIR" -maxdepth 2 -name "scaffolds.fasta" | head -10 | while read f; do
    echo "  - $(dirname "$f" | xargs basename)"
done

# Run pipeline
echo ""
echo "========================================="
echo "Running Pipeline"
echo "========================================="
echo ""
echo "Command:"
echo "  nextflow run main.nf -config nextflow.config \\"
echo "    -profile local_highCPU \\"
echo "    -work-dir '$WORK_DIR' \\"
echo "    --outdir '$OUTDIR' \\"
echo "    --scaffolds '$ASSEMBLIES_DIR' \\"
echo "    --blast_query '$QUERY'"
echo ""

nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir "$WORK_DIR" \
    --outdir "$OUTDIR" \
    --scaffolds "$ASSEMBLIES_DIR" \
    --blast_query "$QUERY"

# Check results
echo ""
echo "========================================="
echo "Checking Results"
echo "========================================="

if [ -d "$OUTDIR/blast_hit_sequences" ]; then
    echo "✓ blast_hit_sequences directory created"
    
    # Count processed samples
    fasta_count=$(ls -1 "$OUTDIR/blast_hit_sequences"/*_top_hits.fasta 2>/dev/null | wc -l)
    tsv_count=$(ls -1 "$OUTDIR/blast_hit_sequences"/*_top_hits_info.tsv 2>/dev/null | wc -l)
    
    echo ""
    echo "Processed samples: $fasta_count FASTA files, $tsv_count TSV files"
    
    # Show samples with hits
    echo ""
    echo "Samples with BLAST hits (non-empty FASTAs):"
    for f in "$OUTDIR/blast_hit_sequences"/*_top_hits.fasta; do
        if [ -s "$f" ]; then
            size=$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")
            sample=$(basename "$f" _top_hits.fasta)
            echo "  - $sample: $size bytes"
        fi
    done | head -20
    
    # Check combined outputs
    echo ""
    echo "Combined outputs:"
    if [ -f "$OUTDIR/all_blast_hits.fasta" ]; then
        size=$(stat -f%z "$OUTDIR/all_blast_hits.fasta" 2>/dev/null || stat -c%s "$OUTDIR/all_blast_hits.fasta")
        seq_count=$(grep -c "^>" "$OUTDIR/all_blast_hits.fasta" || echo "0")
        echo "  ✓ all_blast_hits.fasta: $size bytes, $seq_count sequences"
        
        echo ""
        echo "  First 5 sequence headers:"
        grep "^>" "$OUTDIR/all_blast_hits.fasta" | head -5 | sed 's/^/    /'
    fi
    
    if [ -f "$OUTDIR/all_blast_hits_summary.tsv" ]; then
        lines=$(wc -l < "$OUTDIR/all_blast_hits_summary.tsv")
        data_lines=$((lines - 1))
        echo ""
        echo "  ✓ all_blast_hits_summary.tsv: $data_lines data rows"
        
        echo ""
        echo "  Summary table (first 10 rows):"
        head -11 "$OUTDIR/all_blast_hits_summary.tsv" | column -t -s $'\t' | sed 's/^/    /'
    fi
    
    # Check assembly summary
    echo ""
    echo "Assembly QC summary:"
    if [ -f "$OUTDIR/assembly_summary_metrics.tsv" ]; then
        lines=$(wc -l < "$OUTDIR/assembly_summary_metrics.tsv")
        echo "  ✓ assembly_summary_metrics.tsv: $((lines - 1)) assemblies"
    fi
    
else
    echo "✗ blast_hit_sequences directory NOT found"
    exit 1
fi

echo ""
echo "========================================="
echo "Test Complete!"
echo "========================================="
echo ""
echo "Output directory: $OUTDIR"
echo "Combined FASTA: $OUTDIR/all_blast_hits.fasta"
echo "Combined TSV: $OUTDIR/all_blast_hits_summary.tsv"
