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
include { GTHsearch; ExtractGTHCds; CombineGTHHits } from '../modules/gth_search'
include { AlignMAFFT } from '../modules/align_mafft'
include { AlignMAFFT as AlignGTHCdsMAFFT } from '../modules/align_mafft'
include { AlignMAFFT as AlignProteinMAFFT } from '../modules/align_mafft'
include { MergeBlastSummaries } from '../modules/summarize_blast'
include { PipelineStatistics } from '../modules/pipeline_statistics'

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

    // When GTH is requested with a nucleotide query, verify the query is CDS
    // (no stop codons).  If stop codons are found the query is likely genomic
    // DNA with introns, which GTH cannot handle correctly as an EST input.
    if (params.search_tool == 'gth' && params.blast_query) {
        def qfile = file(params.blast_query)
        if (qfile.exists()) {
            def buf = new StringBuilder()
            def inSeq = false
            qfile.eachLine { line ->
                if (line.startsWith('>')) { inSeq = true }
                else if (inSeq) { buf.append(line.trim().toUpperCase()) }
            }
            def seq = buf.toString()
            if (seq) {
                def nucCount = seq.findAll { c -> 'ATGCUN'.contains(c as String) }.size()
                def isNucl = (nucCount * 100 / seq.size()) >= 95
                if (isNucl) {
                    def stopCodons = ['TAA', 'TAG', 'TGA'] as Set
                    def hasStop = false
                    // Only check frame 0: a CDS must start at position 0.
                    // Checking other frames would produce false positives on any valid CDS.
                    for (int i = 0; i <= seq.size() - 3 && !hasStop; i += 3) {
                        // ignore stop codon at the terminal position
                        if (stopCodons.contains(seq.substring(i, i + 3)) && (i + 3) < seq.size()) hasStop = true
                    }
                    if (hasStop) {
                        exit 1, """ERROR: --search_tool gth requires a CDS or protein query, but the nucleotide query contains stop codons in at least one reading frame.
  This suggests the query is genomic DNA with introns, which GenomeThreader cannot use as an EST/cDNA input.
  Please provide either:
    - A mature CDS sequence (spliced, no introns, no in-frame stop codons)
    - A protein sequence (amino acid FASTA)
  To use genomic DNA as a locus anchor, run without --search_tool gth (BLAST only)."""
                    }
                }
            }
        }
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
                flat: false,
                size: 2
            )
            local_pe_formatted = read_pairs_local_ch.map { sample_id, reads_list ->
                // Get basenames without extensions
                def name1 = reads_list[0].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                def name2 = reads_list[1].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')

                // Find longest common prefix (this is the true sample ID)
                def commonPrefix = ""
                def minLen = Math.min(name1.length(), name2.length())
                for (int i = 0; i < minLen; i++) {
                    if (name1[i] == name2[i]) {
                        commonPrefix += name1[i]
                    } else {
                        break
                    }
                }

                // Remove trailing separators (_, ., -)
                commonPrefix = commonPrefix.replaceAll(/[._-]+$/, '')

                return [commonPrefix, reads_list[0], reads_list[1], 'local']
            }
        } else {
            local_pe_formatted = Channel.empty()
        }

        // Handle SRA paired-end reads
        if (params.SRA_index) {
            sra_pe_formatted = SRAdownloadPE.out.map { sample_id, read1, read2 ->
                [sample_id, read1, read2, 'SRA']
            }
            
            // Handle SRA single-end reads
            sra_se_formatted = SRAdownloadSE.out.map { sample_id, read ->
                [sample_id, read, 'SRA']
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
        if (scaffolds_path.contains('**') || scaffolds_path.contains('*')) {
            // It's a glob pattern - use directly
            local_scaffolds_ch = Channel
                .fromPath(scaffolds_path, checkIfExists: true)
                .map { file ->
                    def base = file.baseName.replaceAll(/\.(fasta|fa|fna)$/, "")
                    // If the filename is a generic name (e.g. "scaffolds"), use
                    // the parent directory name as the sample_id instead
                    def sample_id = (base ==~ /scaffolds?|assembly|genome|contigs?/) ? file.parent.name : base
                    [sample_id, file]
                }
        } else if (file(scaffolds_path).exists() && file(scaffolds_path).isDirectory()) {
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
    // Ensure both channels exist, even if empty
    if (params.reads || params.SRA_index) {
        if (params.scaffolds) {
            all_scaffolds_ch = spades_scaffolds.mix(local_scaffolds_ch)
        } else {
            all_scaffolds_ch = spades_scaffolds
        }
    } else {
        if (params.scaffolds) {
            all_scaffolds_ch = local_scaffolds_ch
        } else {
            all_scaffolds_ch = Channel.empty()
        }
    }

    // ---------------------
    // Rename samples (optional)
    // ---------------------
    if (params.SRR_sample_map) {
        def rename_map = [:]
        file(params.SRR_sample_map).eachLine { line ->
            def parts = line.trim().split(',', 2)
            if (parts.size() == 2) rename_map[parts[0].trim()] = parts[1].trim()
        }
        all_scaffolds_ch = all_scaffolds_ch
            .map { sample_id, scaffolds ->
                [rename_map.getOrDefault(sample_id, sample_id), scaffolds]
            }
    }

    // ---------------------
    // Quality assessment
    // ---------------------
    if (!params.skip_quast) {
        QUAST(all_scaffolds_ch)
        quast_metrics = QUAST.out.metrics.collect()
        SummarizeAssemblyMetrics(quast_metrics)
        assembly_qc_done = SummarizeAssemblyMetrics.out
    } else {
        assembly_qc_done = Channel.value([])
    }

    // ---------------------
    // Sequence search analysis
    // ---------------------
    if (params.blast_query) {
        blast_query_file = file(params.blast_query)

        // Split query file into individual sequences
        SplitQuerySequences(blast_query_file)

        // Read query_info.tsv keyed by filename (exact, reliable key)
        query_info_ch = SplitQuerySequences.out.query_info
            .splitCsv(header: true, sep: '\t')
            .map { row -> [row.Filename, row.Query_ID, row.Sequence_Type] }

        // Join individual query files with their query_id and seq_type
        queries_with_type = SplitQuerySequences.out.individual_queries
            .flatten()
            .map { query_file -> [query_file.name, query_file] }
            .join(query_info_ch)
            .map { filename, query_file, query_id, seq_type ->
                [query_id, query_file, seq_type]
            }

        // [query_id, query_file] — used by both BLAST and GTH paths
        individual_queries = queries_with_type
            .map { query_id, query_file, seq_type -> [query_id, query_file] }

        // Shared input tuple structure for both tools:
        // [sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]
        search_input = all_scaffolds_ch
            .combine(queries_with_type)
            .map { sample_id, scaffolds_fasta, query_id, query_fasta, seq_type ->
                [sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]
            }

        // ------------------------------------------------------------------
        // BLAST search (always runs)
        // tblastn (protein query) or blastn (nucleotide query) against
        // assembled scaffolds, then extract padded hit regions.
        // ------------------------------------------------------------------
        BLASTsearch(search_input)

        // Join BLAST results back with scaffolds for region extraction
        // Format: [query_id, sample_id, blast_results, scaffolds_fasta]
        blast_remapped = BLASTsearch.out.results
            .map { query_id, sample_id, blast_results ->
                [sample_id, [query_id, blast_results]]
            }

        scaffolds_for_join = all_scaffolds_ch
            .map { sample_id, scaffolds -> [sample_id, scaffolds] }

        blast_with_scaffolds = blast_remapped
            .combine(scaffolds_for_join, by: 0)
            .map { sample_id, blast_info, scaffolds_fasta ->
                def query_id      = blast_info[0]
                def blast_results = blast_info[1]
                [query_id, sample_id, blast_results, scaffolds_fasta]
            }

        ExtractBlastHits(blast_with_scaffolds)

        // Group per-sample outputs by query_id
        sequences_by_query = ExtractBlastHits.out.sequences
            .map { query_id, sample_id, fasta -> [query_id, fasta] }
            .groupTuple()

        hit_info_by_query = ExtractBlastHits.out.hit_info
            .groupTuple()

        combined_by_query = individual_queries
            .join(sequences_by_query)
            .join(hit_info_by_query)

        CombineBlastHits(combined_by_query)

        // Merge per-query summary TSVs into a single wide-format table
        all_summaries = CombineBlastHits.out.summary
            .map { query_id, summary_file -> summary_file }
            .collect()

        MergeBlastSummaries(all_summaries)

        // Align BLAST hit sequences across samples
        AlignMAFFT(CombineBlastHits.out.combined_sequences)

        // ------------------------------------------------------------------
        // GenomeThreader path  (--search_tool gth)
        // Runs on BLAST-extracted regions (not the full assembly) for speed.
        // Performs protein-to-genome or cDNA-to-genome spliced alignment to
        // resolve exact exon boundaries. Outputs nucleotide CDS and protein.
        // ------------------------------------------------------------------
        if (params.search_tool == 'gth') {

            // Build GTH input from BLAST-extracted hit FASTAs (not full scaffolds):
            // [sample_id, hit_fasta, query_id, query_fasta, seq_type]
            gth_input = ExtractBlastHits.out.sequences
                .combine(queries_with_type, by: 0)
                .map { query_id, sample_id, hit_fasta, query_fasta, seq_type ->
                    [sample_id, hit_fasta, query_id, query_fasta, seq_type]
                }

            GTHsearch(gth_input)

            // Extract CDS + protein from the raw GFF3 using python3 (spades container)
            ExtractGTHCds(GTHsearch.out.gff3)

            // Group per-sample outputs by query_id
            gth_cds_by_query = ExtractGTHCds.out.cds_sequences
                .map { query_id, sample_id, cds_fasta -> [query_id, cds_fasta] }
                .groupTuple()

            gth_prot_by_query = ExtractGTHCds.out.protein_sequences
                .map { query_id, sample_id, prot_fasta -> [query_id, prot_fasta] }
                .groupTuple()

            gth_summary_by_query = ExtractGTHCds.out.summary
                .groupTuple()

            // Join query file + all per-sample outputs for each query_id
            gth_combined_input = individual_queries
                .join(gth_cds_by_query)
                .join(gth_prot_by_query)
                .join(gth_summary_by_query)

            CombineGTHHits(gth_combined_input)

            // Align nucleotide CDS sequences across samples
            AlignGTHCdsMAFFT(CombineGTHHits.out.combined_cds)

            // Align protein sequences across samples
            AlignProteinMAFFT(CombineGTHHits.out.combined_protein)

            gth_done = CombineGTHHits.out.combined_cds
                .map { query_id, fasta -> fasta }
            all_done = assembly_qc_done
                .mix(MergeBlastSummaries.out)
                .mix(gth_done)
                .collect()

        } else {
            all_done = assembly_qc_done
                .mix(MergeBlastSummaries.out)
                .collect()
        }

    } else {
        all_done = assembly_qc_done.collect()
    }

    // ---------------------
    // Generate pipeline execution statistics
    // ---------------------
    PipelineStatistics(all_done, workflow.launchDir.resolve(params.outdir).resolve("pipeline_trace.txt"))
}
