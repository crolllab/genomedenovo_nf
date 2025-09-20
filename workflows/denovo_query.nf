nextflow.enable.dsl=2

// ---------------------
// Module includes
// ---------------------
include { SRAresolve } from '../modules/resolve_SRA'
include { SRAdownloadPE } from '../modules/download_SRA'
include { SRAdownloadSE } from '../modules/download_SRA'
include { SPAdes } from '../modules/assemble_spades'

// ---------------------
// Main workflow
// ---------------------
workflow variant_calling {

    // ---------------------
    // Input checks
    // ---------------------
    
    if (!params.reads && !params.SRA_index) {
       exit 1, "ERROR: No input reads provided. Supply local FASTQ files and/or SRA run accessions."
    }
    
    // ---------------------
    // SRA metadata + downloads
    // ---------------------
    if (params.SRA_index) {

    // Input SRA index file
    // Step 0: Make a channel from the SRA index file
    sra_file_ch = Channel.fromPath(params.SRA_index)
    
    // First run SRAresolve
    SRAresolve(sra_file_ch)

    // Then use its outputs to create channels for download processes
    pe_ids = SRAresolve.out.pe_file
        .splitText()
        .map { it.trim() }

    se_ids = SRAresolve.out.se_file
        .splitText()
        .map { it.trim() }

    // Now run downloads - these will wait for SRAresolve to complete
    SRAdownloadPE(pe_ids)
    SRAdownloadSE(se_ids)

    }

    // ---------------------
    // Local FASTQ reads
    // ---------------------
if (params.reads) {
    read_pairs_local_ch = Channel.fromFilePairs(
        params.reads,
        checkIfExists: true,
        flat: false // keeps paired R1/R2 in a tuple
    )
    local_pe_formatted = read_pairs_local_ch.map { sample_id, reads_list ->
        [sample_id, reads_list[0], reads_list[1]]
    }
} else {
    local_pe_formatted = Channel.empty()
}

// Only create SRA channel if SRA processing was enabled
if (params.SRA_index) {
    sra_pe_formatted = SRAdownloadPE.out.map { sample_id, read1, read2 ->
        [sample_id, read1, read2]
    }
} else {
    sra_pe_formatted = Channel.empty()
}

combined_pe_ch = sra_pe_formatted.mix(local_pe_formatted)

// ---------------------
// Denovo assembly
// ---------------------
// Connect channels to trimming processes
//trimSequencesPE(combined_pe_ch)

// Handle SE trimming conditionally
if (params.SRA_index) {
//    trimSequencesSE(SRAdownloadSE.out)
} else {
//    trimSequencesSE(Channel.empty())
}
    // Access trimmed outputs
//    pe_trimmed = trimSequencesPE.out.reads
//    se_trimmed = trimSequencesSE.out.reads
//    pe_reports = trimSequencesPE.out.report
//    se_reports = trimSequencesSE.out.report

    // Collect the JSON reports (second output channel, `emit: report`)
//    se_reports_ch = se_reports.collect()
//    pe_reports_ch = pe_reports.collect()


}
