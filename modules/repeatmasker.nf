process RepeatMask {
    tag "RepeatMasker"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/12_repeatmasker", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds)
    path repeat_library

    output:
    tuple val(sample_id), path("*_repeatmasker.gff"),    emit: gff
    tuple val(sample_id), path("*_repeatmasker.masked"), emit: masked
    tuple val(sample_id), path("*_repeatmasker.tbl"),    emit: summary

    script:
    // Use a filesystem-safe ID for file names; original name is preserved as val(sample_id)
    def safe_id = sample_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    set -euo pipefail

    # Copy scaffold to a predictable filename so RepeatMasker output names are known
    cp ${scaffolds} input_scaffolds.fasta

    RepeatMasker -lib "$repeat_library" -pa $task.cpus -gff input_scaffolds.fasta

    mv input_scaffolds.fasta.out.gff "${safe_id}_repeatmasker.gff"
    mv input_scaffolds.fasta.masked  "${safe_id}_repeatmasker.masked"
    mv input_scaffolds.fasta.tbl     "${safe_id}_repeatmasker.tbl"
    """
}
