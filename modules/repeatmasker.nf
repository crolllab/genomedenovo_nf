process RepeatMasker {
    tag "Repeat masking assembly"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/12_repeatmasker", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)
    path repeat_library

    output:
    tuple val(sample_id), path("*_repeatmasker.gff"),    emit: gff
    tuple val(sample_id), path("*_repeatmasker.masked"), emit: masked
    tuple val(sample_id), path("*_repeatmasker.tbl"),    emit: summary
    // The .out table carries the TE class/family (the GFF does not); used by SummarizeHitTEs
    tuple val(sample_id), path("*_repeatmasker.out"),    emit: rm_out

    script:
    // Use a filesystem-safe ID for file names; original name is preserved as val(sample_id)
    def safe_id = sample_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    // -pa sets parallel batches, and each RMBlast batch runs 4 threads
    def batches = Math.max(1, task.cpus.intdiv(4))
    """
    set -euo pipefail

    # Copy scaffold to a predictable filename so RepeatMasker output names are known
    cp ${scaffolds} input_scaffolds.fasta

    RepeatMasker -lib "$repeat_library" -pa ${batches} -gff input_scaffolds.fasta

    mv input_scaffolds.fasta.out.gff "${safe_id}_repeatmasker.gff"
    mv input_scaffolds.fasta.masked  "${safe_id}_repeatmasker.masked"
    mv input_scaffolds.fasta.tbl     "${safe_id}_repeatmasker.tbl"
    mv input_scaffolds.fasta.out     "${safe_id}_repeatmasker.out"
    """
}
