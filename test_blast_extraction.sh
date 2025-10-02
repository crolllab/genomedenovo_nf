#!/bin/bash
# Quick test script to verify the extract_blast_hits.nf fix

echo "========================================="
echo "Testing extract_blast_hits.nf Fix"
echo "========================================="

# Set up test parameters
OUTDIR="./test_blast_extraction"
WORK_DIR="/scratch/daniel_tmp_test"
SCAFFOLDS="./IPO323_sub.fasta"
QUERY="./wagA_gene.fa"

# Clean up previous test
echo "Cleaning up previous test outputs..."
rm -rf "$OUTDIR" "$WORK_DIR"

# Activate environment
echo "Activating micromamba environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf_gp_env

# Run pipeline on pre-assembled scaffolds
echo ""
echo "Running pipeline with scaffolds: $SCAFFOLDS"
echo "BLAST query: $QUERY"
echo ""

nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir "$WORK_DIR" \
    --outdir "$OUTDIR" \
    --scaffolds "$SCAFFOLDS" \
    --blast_query "$QUERY"

# Check results
echo ""
echo "========================================="
echo "Checking Results"
echo "========================================="

if [ -d "$OUTDIR/blast_hit_sequences" ]; then
    echo "✓ blast_hit_sequences directory created"
    
    # Check FASTA files
    echo ""
    echo "FASTA files:"
    for f in "$OUTDIR/blast_hit_sequences"/*.fasta; do
        if [ -f "$f" ]; then
            size=$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")
            lines=$(wc -l < "$f")
            headers=$(grep -c "^>" "$f" || echo "0")
            echo "  - $(basename "$f"): $size bytes, $lines lines, $headers sequences"
        fi
    done
    
    # Check TSV files
    echo ""
    echo "TSV files:"
    for f in "$OUTDIR/blast_hit_sequences"/*.tsv; do
        if [ -f "$f" ]; then
            lines=$(wc -l < "$f")
            data_lines=$((lines - 1))
            echo "  - $(basename "$f"): $lines total lines, $data_lines data rows"
        fi
    done
    
    # Show sample of extracted sequences
    echo ""
    echo "Sample of extracted sequences:"
    head -20 "$OUTDIR/blast_hit_sequences"/*_top_hits.fasta
    
    # Check combined outputs
    echo ""
    echo "Combined outputs:"
    if [ -f "$OUTDIR/all_blast_hits.fasta" ]; then
        size=$(stat -f%z "$OUTDIR/all_blast_hits.fasta" 2>/dev/null || stat -c%s "$OUTDIR/all_blast_hits.fasta")
        headers=$(grep -c "^>" "$OUTDIR/all_blast_hits.fasta" || echo "0")
        echo "  - all_blast_hits.fasta: $size bytes, $headers sequences"
    fi
    
    if [ -f "$OUTDIR/all_blast_hits_summary.tsv" ]; then
        lines=$(wc -l < "$OUTDIR/all_blast_hits_summary.tsv")
        echo "  - all_blast_hits_summary.tsv: $lines lines"
        echo ""
        echo "    Contents:"
        cat "$OUTDIR/all_blast_hits_summary.tsv"
    fi
    
else
    echo "✗ blast_hit_sequences directory NOT found"
    exit 1
fi

echo ""
echo "========================================="
echo "Test Complete"
echo "========================================="
