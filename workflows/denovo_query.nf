nextflow.enable.dsl=2

// ---------------------
// Module includes
// ---------------------
include { SRAresolve } from '../modules/resolve_SRA'
include { SRAdownloadPE; SRAdownloadSE } from '../modules/download_SRA'
include { SPAdes; SPAdesSE } from '../modules/assemble_spades'
include { QUAST; SummarizeAssemblyMetrics } from '../modules/assembly_qc'
include { SplitQuerySequences; BLASTsearch } from '../modules/blast_search'
include { ExtractBlastHits; CombineBlastHits } from '../modules/extract_blast_hits'
include { AlignMAFFT } from '../modules/align_mafft'
include { MergeBlastSummaries } from '../modules/summarize_blast'

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
        
        // First run SRAresolve to get URLs
        SRAresolve(sra_file_ch)
        
        // Parse the URL file to create download tuples
        pe_download_info = SRAresolve.out.url_file
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.Layout == 'PE' }
            .map { row -> [row.SRR_ID, row.URL_1, row.URL_2, row.Source] }
        
        se_download_info = SRAresolve.out.url_file
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.Layout == 'SE' }
            .map { row -> [row.SRR_ID, row.URL_1, row.Source] }
        
        // Now run hybrid downloads (ENA first, fallback to NCBI)
        SRAdownloadPE(pe_download_info)
        SRAdownloadSE(se_download_info)
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

        // Handle SRA paired-end reads
        if (params.SRA_index) {
            sra_pe_formatted = SRAdownloadPE.out.map { sample_id, read1, read2 ->
                [sample_id, read1, read2]
            }
            
            // Handle SRA single-end reads
            sra_se_formatted = SRAdownloadSE.out.map { sample_id, read ->
                [sample_id, read]
            }
        } else {
            sra_pe_formatted = Channel.empty()
            sra_se_formatted = Channel.empty()
        }

        // Combine PE reads from different sources
        combined_pe_ch = sra_pe_formatted.mix(local_pe_formatted)

        // Run SPAdes assembly for paired-end
        SPAdes(combined_pe_ch)
        
        // Run SPAdes assembly for single-end
        SPAdesSE(sra_se_formatted)
        
        // Combine scaffolds from both PE and SE assemblies
        spades_scaffolds = SPAdes.out.scaffolds.mix(SPAdesSE.out.scaffolds)
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
    if (params.blast_query) {
        blast_query_file = file(params.blast_query)
        
        // Split query file into individual sequences
        SplitQuerySequences(blast_query_file)
        
        // Read query_info.tsv to get sequence types for each query
        query_info_ch = SplitQuerySequences.out.query_info
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                // Map: [query_id, seq_type]
                [row.Query_ID, row.Sequence_Type]
            }
        
        // Create channel of individual query files with query IDs
        individual_queries = SplitQuerySequences.out.individual_queries
            .flatten()
            .map { query_file ->
                // Extract query ID from filename: query_1_seqname.fasta -> seqname
                def filename = query_file.name
                def query_id = filename.replaceAll(/^query_\d+_/, '').replaceAll(/\.fasta$/, '')
                [query_id, query_file]
            }
        
        // Join query files with their sequence types
        queries_with_type = individual_queries
            .join(query_info_ch)
            .map { query_id, query_file, seq_type ->
                [query_id, query_file, seq_type]
            }
        
        // Combine scaffolds with each query for BLAST
        // Creates: [sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]
        blast_input = all_scaffolds_ch
            .combine(queries_with_type)
            .map { sample_id, scaffolds_fasta, query_id, query_fasta, seq_type ->
                [sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]
            }
        
        // Run BLAST on all scaffolds vs each query
        BLASTsearch(blast_input)
        
        // Join BLAST results with scaffold files for extraction
        // Remap BLAST results to [sample_id, [query_id, blast_results]]
        // and join with [sample_id, scaffolds_fasta]
        // Format: [query_id, sample_id, blast_results, scaffolds_fasta]
        blast_remapped = BLASTsearch.out.results
            .map { query_id, sample_id, blast_results -> 
                [sample_id, [query_id, blast_results]]
            }
        
        scaffolds_for_join = all_scaffolds_ch.map { sample_id, scaffolds -> 
            [sample_id, scaffolds]
        }
        
        blast_with_scaffolds = blast_remapped
            .combine(scaffolds_for_join, by: 0)
            .map { sample_id, blast_info, scaffolds_fasta ->
                def query_id = blast_info[0]
                def blast_results = blast_info[1]
                [query_id, sample_id, blast_results, scaffolds_fasta]
            }
        
        // Extract hit regions from assemblies
        ExtractBlastHits(blast_with_scaffolds)
        
        // Group by query_id for combining
        sequences_by_query = ExtractBlastHits.out.sequences
            .map { query_id, sample_id, fasta -> [query_id, fasta] }
            .groupTuple()
        
        hit_info_by_query = ExtractBlastHits.out.hit_info
            .groupTuple()
        
        // Join with query sequences to include in output
        queries_for_combine = individual_queries
            .map { query_id, query_file -> [query_id, query_file] }
        
        // Combine sequences and hit info for each query
        combined_by_query = queries_for_combine
            .join(sequences_by_query)
            .join(hit_info_by_query)
        
        // Combine all extracted sequences for each query (with query sequence prepended)
        CombineBlastHits(combined_by_query)
        
        // ---------------------
        // Merge BLAST summaries
        // ---------------------
        // Collect all summary TSV files and merge into a single wide-format table
        all_summaries = CombineBlastHits.out.summary
            .map { query_id, summary_file -> summary_file }
            .collect()
        
        MergeBlastSummaries(all_summaries)
        
        // ---------------------
        // Multiple sequence alignment
        // ---------------------
        // Run MAFFT alignment on each query's combined BLAST hit sequences
        AlignMAFFT(CombineBlastHits.out.combined_sequences)
    }
}
