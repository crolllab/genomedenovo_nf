nextflow.enable.dsl=2

// ---------------------
// Module includes
// ---------------------
include { SRAresolve } from '../modules/resolve_SRA'
include { SRAdownloadPE; SRAdownloadSE } from '../modules/download_SRA'
include { SPAdes } from '../modules/assemble_spades'
include { QUAST; SummarizeAssemblyMetrics } from '../modules/assembly_qc'
include { BLASTn } from '../modules/blast_search'
include { ExtractBlastHits; CombineBlastHits } from '../modules/extract_blast_hits'

// ---------------------
// Main workflow
// ---------------------
workflow denovo_query {

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
    SPAdes(combined_pe_ch)
    // Run QUAST on scaffolds
    QUAST(SPAdes.out.scaffolds)

    // Collect all QUAST metrics and summarize
    quast_metrics = QUAST.out.metrics.collect()
    SummarizeAssemblyMetrics(quast_metrics)

    blast_query_file = file(params.blast_query)
    BLASTn(SPAdes.out.scaffolds, blast_query_file)

    // Extract hit regions from assemblies
    ExtractBlastHits(
        BLASTn.out.results,
        SPAdes.out.scaffolds
    )
    
    // Combine all extracted sequences
    all_hit_sequences = ExtractBlastHits.out.sequences.map { it[1] }.collect()
    all_hit_info = ExtractBlastHits.out.hit_info.collect()
    
    CombineBlastHits(all_hit_sequences, all_hit_info)

}
