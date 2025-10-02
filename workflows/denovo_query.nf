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
    if (!params.reads && !params.SRA_index && !params.scaffolds) {
        exit 1, "ERROR: No input provided. Supply local FASTQ files, SRA run accessions, or local scaffold files."
    }
    
    if (params.blast_query && !file(params.blast_query).exists()) {
        exit 1, "ERROR: BLAST query file does not exist: ${params.blast_query}"
    }

    // ---------------------
    // SRA metadata + downloads
    // ---------------------
    if (params.SRA_index) {
        // Input SRA index file
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
        
        // Now run downloads
        SRAdownloadPE(pe_ids)
        SRAdownloadSE(se_ids)
    }

    // ---------------------
    // Local FASTQ reads and assembly
    // ---------------------
    if (params.reads || params.SRA_index) {
        // Handle local FASTQ reads
        if (params.reads) {
            read_pairs_local_ch = Channel.fromFilePairs(
                params.reads,
                checkIfExists: true,
                flat: false
            )
            local_pe_formatted = read_pairs_local_ch.map { sample_id, reads_list ->
                [sample_id, reads_list[0], reads_list[1]]
            }
        } else {
            local_pe_formatted = Channel.empty()
        }

        // Handle SRA reads
        if (params.SRA_index) {
            sra_pe_formatted = SRAdownloadPE.out.map { sample_id, read1, read2 ->
                [sample_id, read1, read2]
            }
        } else {
            sra_pe_formatted = Channel.empty()
        }

        combined_pe_ch = sra_pe_formatted.mix(local_pe_formatted)

        // Run SPAdes assembly
        SPAdes(combined_pe_ch)
        spades_scaffolds = SPAdes.out.scaffolds
    } else {
        spades_scaffolds = Channel.empty()
    }

    // ---------------------
    // Local scaffold files
    // ---------------------
    if (params.scaffolds) {
        scaffolds_path = params.scaffolds
        
        // Check if it's a directory containing subdirectories with scaffolds.fasta
        if (file(scaffolds_path).isDirectory()) {
            // Pattern to find scaffolds.fasta in subdirectories: parent_dir/*/scaffolds.fasta
            scaffolds_pattern = "${scaffolds_path}/*/scaffolds.{fasta,fa,fna}"
            
            local_scaffolds_ch = Channel
                .fromPath(scaffolds_pattern, checkIfExists: true)
                .map { file ->
                    // Use the parent directory name as sample_id
                    def sample_id = file.parent.name
                    [sample_id, file]
                }
        } else {
            // If it's a file pattern or specific file(s)
            local_scaffolds_ch = Channel
                .fromPath(scaffolds_path, checkIfExists: true)
                .map { file ->
                    def sample_id = file.baseName.replaceAll(/\.(fasta|fa|fna)$/, "")
                    [sample_id, file]
                }
        }
    } else {
        local_scaffolds_ch = Channel.empty()
    }

    // ---------------------
    // Combine all scaffolds for downstream analysis
    // ---------------------
    all_scaffolds_ch = spades_scaffolds.mix(local_scaffolds_ch)

    // ---------------------
    // Quality assessment
    // ---------------------
    // Run QUAST on all scaffolds (only if we have any)
    QUAST(all_scaffolds_ch)
    
    // Collect all QUAST metrics and summarize
    quast_metrics = QUAST.out.metrics.collect()
    SummarizeAssemblyMetrics(quast_metrics)

    // ---------------------
    // BLAST analysis
    // ---------------------
// BLAST analysis
if (params.blast_query) {
    blast_query_file = file(params.blast_query)
    
    // Run BLAST on all scaffolds
    BLASTn(all_scaffolds_ch, blast_query_file)
    
     // Join BLAST results with scaffold files for extraction
    blast_with_scaffolds = BLASTn.out.results.join(all_scaffolds_ch)
    
    // Extract hit regions from assemblies
    ExtractBlastHits(blast_with_scaffolds)
    
    // Combine all extracted sequences
    all_hit_sequences = ExtractBlastHits.out.sequences.map { it[1] }.collect()
    all_hit_info = ExtractBlastHits.out.hit_info.collect()
    
    CombineBlastHits(all_hit_sequences, all_hit_info)
}
}
