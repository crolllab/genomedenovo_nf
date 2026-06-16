// Align single-copy BUSCO ortholog sequences across samples using MAFFT.
// Input: one [gene_id, [faa_files_from_all_samples]] tuple per gene.
// Output: a single aligned FASTA per gene for downstream supermatrix construction.
process AlignBuscoGene {
    tag "MAFFT alignment"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

    input:
    tuple val(gene_id), path(faa_files, stageAs: 'sample_??.faa')

    output:
    path "${gene_id}_aligned.fasta", emit: alignment

    script:
    """
    set -euo pipefail

    # Concatenate all per-sample sequences for this gene into one multi-FASTA
    # (stageAs gives each file a unique name to avoid staging collisions)
    cat sample_*.faa > ${gene_id}_combined.fasta

    # Align with MAFFT auto strategy (--amino for protein sequences)
    mafft --amino --auto --thread ${task.cpus} ${gene_id}_combined.fasta > ${gene_id}_aligned.fasta
    """
}
