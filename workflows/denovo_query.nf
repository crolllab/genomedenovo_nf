// ---------------------
// Module includes
// ---------------------
include { SRAresolve } from '../modules/resolve_SRA'
include { SRAdownloadPE; SRAdownloadSE } from '../modules/download_SRA'
include { trimSequencesPE; trimSequencesSE } from '../modules/fastp_trimming'
include { SPAdes; SPAdesSE } from '../modules/assemble_spades'
include { QUAST; SummarizeAssemblyMetrics } from '../modules/assembly_qc'
include { DownloadBUSCODatabase; BUSCO; SummarizeBuscoMetrics } from '../modules/busco_qc'
include { SplitQuerySequences; BLASTsearch } from '../modules/blast_search'
include { ExtractBlastHits; CombineBlastHits } from '../modules/extract_blast_hits'
include { ExonerateSearch; ExtractExonerateCds; CombineExonerateHits } from '../modules/exonerate_search'
include { AssessCdsIntegrity } from '../modules/cds_integrity'
include { AlignMAFFT } from '../modules/align_mafft'
include { AlignMAFFT as AlignExonerateCdsMAFFT } from '../modules/align_mafft'
include { AlignMAFFT as AlignProteinMAFFT } from '../modules/align_mafft'
include { MergeBlastSummaries } from '../modules/summarize_blast'
include { MergeExonerateSummaries } from '../modules/summarize_exonerate'
include { SummarizeIsoforms } from '../modules/summarize_isoforms'
include { SummarizeFastpStats } from '../modules/summarize_fastp'
include { SummarizeBwaStats } from '../modules/summarize_bwa'
include { bwaIndex } from '../modules/bwa_index'
include { bwaMap } from '../modules/bwa_mapping'
include { samtoolsSort } from '../modules/samtools_sort'
include { RepeatMasker } from '../modules/repeatmasker'
include { AlignBuscoGene } from '../modules/align_busco_genes'
include { BuildSupermatrix } from '../modules/build_supermatrix'
include { IQtreePhylogeny } from '../modules/iqtree_phylogeny'
include { RaxmlNgPhylogeny } from '../modules/raxml_ng_phylogeny'
include { IQtreeAlignmentPhylogeny as IQtreeBlastPhylogeny } from '../modules/iqtree_phylogeny'
include { IQtreeAlignmentPhylogeny as IQtreeCdsPhylogeny } from '../modules/iqtree_phylogeny'
include { IQtreeAlignmentPhylogeny as IQtreeProteinPhylogeny } from '../modules/iqtree_phylogeny'
include { IQtreeAlignmentPhylogeny as IQtreeCdsIsoformPhylogeny } from '../modules/iqtree_phylogeny'
include { IQtreeAlignmentPhylogeny as IQtreeProteinIsoformPhylogeny } from '../modules/iqtree_phylogeny'
include { PlotIsoformTree } from '../modules/plot_isoform_tree'
include { SummarizeHitTEs } from '../modules/summarize_hit_tes'
include { BLASTparalogs; ParalogCheck } from '../modules/paralog_check'
include { BlastLocusOnReference; ComputeLocusCoverage; MergeLocusCoverage } from '../modules/locus_coverage'
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

    if (params.reference_genome && !file(params.reference_genome).exists()) {
        exit 1, "ERROR: Reference genome file does not exist: ${params.reference_genome}"
    }

    if (params.repeat_library && !file(params.repeat_library).exists()) {
        exit 1, "ERROR: Repeat library file does not exist: ${params.repeat_library}"
    }

    if (params.phylo_tree && !params.busco_lineage) {
        exit 1, "ERROR: --phylo_tree requires --busco_lineage to be set."
    }

    // When a pre-downloaded database is supplied, validate the lineage exists inside it
    // before any work is submitted — this catches typos (e.g. odb12 vs odb10) immediately.
    if (params.busco_lineage && params.busco_db_path) {
        def lineage_dir = file("${params.busco_db_path}/lineages/${params.busco_lineage}")
        if (!lineage_dir.exists()) {
            def available = file("${params.busco_db_path}/lineages").exists()
                ? file("${params.busco_db_path}/lineages").list()?.sort()?.join(', ') ?: 'none found'
                : 'lineages/ directory not found'
            exit 1, """ERROR: BUSCO lineage '${params.busco_lineage}' not found in '${params.busco_db_path}/lineages/'.
  Available lineages: ${available}
  Browse all lineages at: https://busco-data.ezlab.org/v5/data/lineages/"""
        }
    }

    if (params.phylo_tree && !(params.busco_tree in ['raxml-ng', 'iqtree'])) {
        exit 1, "ERROR: --busco_tree must be 'raxml-ng' (default) or 'iqtree'. Got: '${params.busco_tree}'."
    }

    if (params.phylo_tree && params.busco_subsample != 'all') {
        if (!(params.busco_subsample instanceof Number) && !(params.busco_subsample ==~ /^[1-9][0-9]*/)) {
            exit 1, "ERROR: --busco_subsample must be a positive integer or 'all'. Got: '${params.busco_subsample}'."
        }
        if ((params.busco_subsample as Integer) < 1) {
            exit 1, "ERROR: --busco_subsample must be >= 1. Got: '${params.busco_subsample}'."
        }
    }

    if (!(params.search_tool in ['blast', 'exonerate'])) {
        exit 1, "ERROR: Unsupported --search_tool '${params.search_tool}'. Allowed values: blast, exonerate."
    }

    if (params.locus_coverage && !params.reference_genome) {
        exit 1, "ERROR: --locus_coverage requires --reference_genome (for BWA mapping and BLAST-based locus search)."
    }

    if (params.locus_coverage && !params.blast_query) {
        exit 1, "ERROR: --locus_coverage requires --blast_query (the query is BLASTed against the reference to identify the locus)."
    }

    if (params.locus_coverage && !params.reads && !params.SRA_index) {
        exit 1, "ERROR: --locus_coverage requires read input (--reads or --SRA_index) to produce BAM files for depth analysis."
    }

    def min_align_prop = params.min_align_prop.toString()
    if (!(min_align_prop ==~ /^[0-9]*\.?[0-9]+$/) || (min_align_prop as BigDecimal) > 1) {
        exit 1, "ERROR: --min_align_prop must be a number between 0 and 1. Got: '${params.min_align_prop}'."
    }

    // --query_paralog: paralogs of the queries, each named '<query>_paralogN'
    // after the first word of a --blast_query header, and of the same type
    if (params.query_paralog) {
        if (!params.blast_query) {
            exit 1, "ERROR: --query_paralog requires --blast_query."
        }
        if (!file(params.query_paralog).exists()) {
            exit 1, "ERROR: Paralog query file does not exist: ${params.query_paralog}"
        }
        // [first word of header, nucleotide or protein] per record
        def fastaRecords = { path ->
            def recs = []
            def name = null
            def buf = new StringBuilder()
            def flush = {
                if (name != null) {
                    def s = buf.toString().toUpperCase()
                    def nuc = s.findAll { c -> 'ATGCUN'.contains(c as String) }.size()
                    recs << [name, (s && nuc * 100 / s.size() >= 95) ? 'nucleotide' : 'protein']
                }
            }
            file(path).eachLine { line ->
                if (line.startsWith('>')) {
                    flush()
                    name = line.substring(1).trim().tokenize(' \t')[0]
                    buf = new StringBuilder()
                } else {
                    buf.append(line.trim())
                }
            }
            flush()
            recs
        }
        def query_types = fastaRecords(params.blast_query).collectEntries { name, type -> [name, type] }
        def paralog_recs = fastaRecords(params.query_paralog)
        def problems = []
        paralog_recs.each { name, type ->
            def query = name.replaceAll(/_paralog[0-9]+$/, '')
            if (query == name) {
                problems << "'${name}': not named '<query>_paralogN'"
            } else if (!query_types.containsKey(query)) {
                problems << "'${name}': no query '${query}' in --blast_query"
            } else if (query_types[query] != type) {
                problems << "'${name}': ${type} sequence, but query '${query}' is ${query_types[query]}"
            }
        }
        def dups = paralog_recs.countBy { name, _type -> name }.findAll { _name, n -> n > 1 }.keySet()
        if (dups) problems << "duplicated names: ${dups.join(', ')}"
        if (!paralog_recs) problems << "no sequences found"
        if (problems) {
            exit 1, """ERROR: Invalid --query_paralog file ${params.query_paralog}:
${problems.collect { p -> '    ' + p }.join('\n')}
  Name each paralog after the first word of its query header plus '_paralog1', '_paralog2', ...
  (e.g. '>GeneA_paralog1' for the query '>GeneA ...'); paralogs must be of the same type as the query."""
        }
        def with_paralogs = paralog_recs.collect { name, _type -> name.replaceAll(/_paralog[0-9]+$/, '') } as Set
        def without = query_types.keySet().findAll { q -> !with_paralogs.contains(q) }
        if (without) {
            log.warn "No paralogs given for query/queries ${without.join(', ')}: no paralog check for them."
        }
    }

    // When Exonerate is requested with a nucleotide query, verify the query is CDS
    // (no stop codons). If stop codons are found the query is likely genomic
    // DNA with introns, which Exonerate cannot handle correctly with coding2genome.
    if (params.search_tool == 'exonerate' && params.blast_query) {
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
                    // Only check frame 0: a CDS must start at position 0.
                    // Checking other frames would produce false positives on any valid CDS.
                    // Exclude the terminal codon (position seq.size()-3 onward)
                    def hasStop = (0..<(seq.size() - 3)).step(3).any { i ->
                        stopCodons.contains(seq.substring(i, i + 3))
                    }
                    if (hasStop) {
                                                exit 1, """ERROR: --search_tool exonerate requires a CDS or protein query, but the nucleotide query contains stop codons in at least one reading frame.
    This suggests the query is genomic DNA with introns, which Exonerate cannot use correctly with the coding2genome model.
  Please provide either:
    - A mature CDS sequence (spliced, no introns, no in-frame stop codons)
    - A protein sequence (amino acid FASTA)
    To use genomic DNA as a locus anchor, run without --search_tool exonerate (BLAST only)."""
                    }
                }
            }
        }
    }

    // ---------------------
    // Sample names
    // ---------------------
    // Inside the pipeline every sample is named by its filesystem-safe name:
    // characters outside [A-Za-z0-9._-] become '_' (e.g. 'Iran_2024_Kch(2)12-1g'
    // -> 'Iran_2024_Kch_2_12-1g'). Sample names end up in file names, FASTA
    // headers and tree tips, and characters such as '(' or ' ' broke steps
    // silently. <outdir>/sample_names.tsv maps every input to its final name.
    def safeName = { name -> name.toString().replaceAll(/[^A-Za-z0-9._-]/, '_') }

    // Optional renaming (--SRR_sample_map), looked up by the original input ID
    def rename_map = [:]
    if (params.SRR_sample_map) {
        file(params.SRR_sample_map).eachLine { line ->
            def parts = line.trim().split(',', 2)
            if (parts.size() == 2) rename_map[parts[0].trim()] = parts[1].trim()
        }
    }
    // Local read samples are already safe-named when they are renamed, so the
    // map is also consulted with the safe form of the original ID
    def renamedTo = { id -> rename_map[id] ?: rename_map.find { k, _v -> safeName(k) == id }?.value ?: id }
    def finalName = { id -> safeName(renamedTo(id)) }

    // Placeholder for fastp summary signal (only populated when reads are trimmed)
    fastp_summary_done = channel.empty()

    // Placeholder for BWA mapping signal (only populated when --reference_genome is set)
    bwa_done = channel.empty()

    // Placeholder for BAM channel passed to locus coverage (set inside the BWA block)
    bam_ch_for_locus = channel.empty()

    // Placeholder for locus coverage signal (only populated when --locus_coverage is set)
    locus_cov_done = channel.empty()

    // Placeholder for RepeatMasker signal (only populated when --repeat_library is set)
    repeat_done = channel.empty()

    // Placeholders for the TEs-near-hits step (--repeat_library with --blast_query)
    repeat_out_all = channel.empty()
    hit_tes_done = channel.empty()

    // Placeholder for BUSCO phylogeny signal (only populated when --phylo_tree is set)
    phylo_done = channel.empty()

    // Placeholder for the paralog check signal (only populated when --query_paralog is set)
    paralog_done = channel.empty()

    // Sample IDs of all inputs (SRA runs, local reads, local scaffolds), known
    // before assembly; used for the sample name collision check
    input_ids_ch = channel.empty()

    // ---------------------
    // SRA metadata + downloads
    // ---------------------
    if (params.SRA_index) {
        // Input SRA index file
        sra_file_ch = channel.fromPath(params.SRA_index)
        
        // First run SRAresolve to get URLs
        SRAresolve(sra_file_ch)
        input_ids_ch = input_ids_ch.mix(
            SRAresolve.out.url_file
                .splitCsv(header: true, sep: '\t')
                .map { row -> row.SRR_ID }
        )

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
            // Support comma-separated list of glob patterns
            def reads_patterns = params.reads.tokenize(';').collect { it.trim() }
            read_pairs_local_ch = channel.fromFilePairs(
                reads_patterns,
                checkIfExists: true,
                flat: false,
                size: 2
            )
            local_pe_raw = read_pairs_local_ch.map { sample_id, reads_list ->
                // Get basenames without extensions
                def name1 = reads_list[0].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                def name2 = reads_list[1].name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')

                // Find longest common prefix (this is the true sample ID)
                def minLen = Math.min(name1.length(), name2.length())
                def commonPrefix = (0..<minLen)
                    .takeWhile { i -> name1[i] == name2[i] }
                    .collect { i -> name1[i] }
                    .join('')

                // Remove trailing separators (_, ., -)
                commonPrefix = commonPrefix.replaceAll(/[._-]+$/, '')

                return [commonPrefix, reads_list[0], reads_list[1], 'local']
            }
            input_ids_ch = input_ids_ch.mix(local_pe_raw.map { row -> row[0] })
            // Trimming and assembly already use the safe name
            local_pe_formatted = local_pe_raw.map { sample_id, r1, r2, source -> [safeName(sample_id), r1, r2, source] }
        } else {
            local_pe_formatted = channel.empty()
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
            sra_pe_formatted = channel.empty()
            sra_se_formatted = channel.empty()
        }

        // Combine PE reads from different sources
        combined_pe_ch = sra_pe_formatted.mix(local_pe_formatted)

        // Trim PE reads with fastp before assembly
        trimSequencesPE(combined_pe_ch)
        trimmed_pe_ch = trimSequencesPE.out.reads.map { sample_id, reads ->
            [sample_id, reads[0], reads[1], 'trimmed']
        }

        // Trim SE reads with fastp before assembly
        trimSequencesSE(sra_se_formatted)
        trimmed_se_ch = trimSequencesSE.out.reads.map { sample_id, read ->
            [sample_id, read, 'trimmed']
        }

        // Collect all fastp JSON reports and generate summary statistics
        fastp_json_ch = trimSequencesPE.out.report
            .mix(trimSequencesSE.out.report)
            .collect()
        fastp_summary_script = Channel.value(file("${projectDir}/bin/summarize_fastp.R"))
        SummarizeFastpStats(fastp_json_ch, fastp_summary_script)
        fastp_summary_done = SummarizeFastpStats.out.tsv

        // ---------------------
        // Optional: BWA-mem2 reference mapping (triggered by --reference_genome)
        // ---------------------
        if (params.reference_genome) {
            ref_ch = Channel.value(file(params.reference_genome))

            // Build BWA-mem2 index once for the reference
            (bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_0123) = bwaIndex(ref_ch)

            // Combine PE and SE trimmed reads into a single channel.
            // SE reads are wrapped in a list so trimmed_reads.size() works
            // correctly (1 vs 2) inside bwaMap for the reads_cmd construction.
            bwa_pe_ch = trimSequencesPE.out.reads
            bwa_se_ch = trimSequencesSE.out.reads.map { id, r -> [id, [r]] }
            bwa_reads_ch = bwa_pe_ch.mix(bwa_se_ch)

            bwaMap(ref_ch, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_0123, bwa_reads_ch)
            samtoolsSort(bwaMap.out.sam)
            bam_ch_for_locus = samtoolsSort.out.bam
            bwa_done = samtoolsSort.out.bam.map { id, bam, bai -> bam }.collect()

            bwa_flagstat_ch = samtoolsSort.out.report.collect()
            bwa_summary_script = Channel.value(file("${projectDir}/bin/summarize_bwa.R"))
            SummarizeBwaStats(bwa_flagstat_ch, bwa_summary_script)
        }

        // Run SPAdes assembly for paired-end
        SPAdes(trimmed_pe_ch)
        
        // Run SPAdes assembly for single-end
        SPAdesSE(trimmed_se_ch)
        
        // Combine scaffolds from both PE and SE assemblies
        spades_scaffolds = SPAdes.out.scaffolds.mix(SPAdesSE.out.scaffolds)
    } else {
        spades_scaffolds = channel.empty()
    }

    // ---------------------
    // Local scaffold files
    // ---------------------
    if (params.scaffolds) {
        // Support comma-separated list of paths/globs (or a single entry)
        def scaffold_entries = params.scaffolds.tokenize(';').collect { it.trim() }

        local_scaffolds_ch = scaffold_entries.inject(channel.empty()) { acc, scaffolds_path ->
            def path_ch
            // Check if it's a glob pattern
            if (scaffolds_path.contains('**') || scaffolds_path.contains('*')) {
                // It's a glob pattern - use directly
                path_ch = channel
                    .fromPath(scaffolds_path, checkIfExists: true)
                    .map { f ->
                        def base = f.baseName.replaceAll(/\.(fasta|fa|fna)$/, "")
                        // If the filename is a generic name (e.g. "scaffolds"), use
                        // the parent directory name as the sample_id instead
                        def sample_id = (base ==~ /scaffolds?|assembly|genome|contigs?/) ? f.parent.name : base
                        [sample_id, f]
                    }
            } else if (file(scaffolds_path).exists() && file(scaffolds_path).isDirectory()) {
                // Pattern to find scaffolds.fasta in subdirectories: parent_dir/*/scaffolds.fasta
                def scaffolds_pattern = "${scaffolds_path}/*/scaffolds.{fasta,fa,fna}"
                path_ch = channel
                    .fromPath(scaffolds_pattern, checkIfExists: true)
                    .map { f ->
                        // Use the parent directory name as sample_id
                        def sample_id = f.parent.name
                        [sample_id, f]
                    }
            } else {
                // If it's a file pattern or specific file(s)
                path_ch = channel
                    .fromPath(scaffolds_path, checkIfExists: true)
                    .map { f ->
                        def sample_id = f.baseName.replaceAll(/\.(fasta|fa|fna)$/, "")
                        [sample_id, f]
                    }
            }
            acc.mix(path_ch)
        }
        input_ids_ch = input_ids_ch.mix(local_scaffolds_ch.map { sample_id, _scaffolds -> sample_id })
    } else {
        local_scaffolds_ch = channel.empty()
    }

    // ---------------------
    // Combine all scaffolds for downstream analysis
    // ---------------------
    // Ensure both channels exist, even if empty
    if (params.reads || params.SRA_index) {
        if (params.scaffolds) {
            // Mix SPAdes and local scaffolds; sample IDs are distinct so no deduplication needed.
            // (Previously used groupTuple() here, which blocked all local scaffolds until
            // every SPAdes job finished — causing QUAST/BLAST to never start on pre-assembled inputs.)
            all_scaffolds_ch = spades_scaffolds.mix(local_scaffolds_ch)
        } else {
            all_scaffolds_ch = spades_scaffolds
        }
    } else {
        if (params.scaffolds) {
            all_scaffolds_ch = local_scaffolds_ch
        } else {
            all_scaffolds_ch = channel.empty()
        }
    }

    // ---------------------
    // Final sample names: --SRR_sample_map renaming, then the safe form
    // ---------------------
    all_scaffolds_ch = all_scaffolds_ch
        .map { sample_id, scaffolds -> [finalName(sample_id), scaffolds] }

    // ---------------------
    // Sample name collision check and name table
    // ---------------------
    // Every sample needs a unique final name. Otherwise outputs overwrite each
    // other and IQ-TREE rejects duplicated names. Runs on the input IDs, so the
    // run stops before assembly (SRA runs are known once SRAresolve has finished).
    input_ids_ch
        .map { input_id -> [finalName(input_id), input_id] }
        .collect(flat: false)
        .map { pairs ->
            def clashes = pairs
                .groupBy { name, _orig -> name }
                .findAll { _name, group -> group.size() > 1 }
            if (clashes) {
                def lines = clashes.collect { name, group ->
                    "    ${name}: " + group.collect { _name, orig ->
                        renamedTo(orig) == orig ? "'${orig}'" : "'${orig}' (renamed to '${renamedTo(orig)}')"
                    }.join(', ')
                }
                error """ERROR: Sample name collision, several inputs end up with the same sample name:
${lines.join('\n')}
  Sample names must be unique, also after characters outside [A-Za-z0-9._-] are replaced by '_'.
  Check --SRR_sample_map for duplicated target names, and inputs (reads, SRA runs, scaffolds)
  given more than once or sharing a name."""
            }
            "input_id\trenamed_to\tsample_id\n" + pairs
                .sort { a, b -> a[0] <=> b[0] }
                .collect { name, orig -> "${orig}\t${renamedTo(orig)}\t${name}\n" }
                .join('')
        }
        .collectFile(name: 'sample_names.tsv', storeDir: params.outdir)

    // ---------------------
    // RepeatMasker TE annotation (optional — triggered by --repeat_library)
    // ---------------------
    if (params.repeat_library) {
        repeat_lib_ch = Channel.value(file(params.repeat_library))
        RepeatMasker(all_scaffolds_ch, repeat_lib_ch)
        repeat_done = RepeatMasker.out.gff.map { id, gff -> gff }.collect()

        // All samples' .out tables as a single list item; wrapped so that
        // `combine` with the per-query hit tables does not flatten it
        repeat_out_all = RepeatMasker.out.rm_out.map { id, rm_out -> rm_out }.collect().map { files -> [files] }
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
        assembly_qc_done = channel.value([])
    }

    // ---------------------
    // BUSCO completeness (optional — only runs when --busco_lineage is set)
    // ---------------------
    if (params.busco_lineage) {
        // Download the lineage database once; if busco_db_path is pre-supplied, use it directly
        if (params.busco_db_path) {
            busco_db_ch = Channel.value(file(params.busco_db_path))
        } else {
            DownloadBUSCODatabase()
            busco_db_ch = DownloadBUSCODatabase.out.db_dir
        }
        BUSCO(all_scaffolds_ch, busco_db_ch)
        SummarizeBuscoMetrics(BUSCO.out.summary.collect())

        // ----- BUSCO gene phylogeny (optional — triggered by --phylo_tree) -----
        if (params.phylo_tree) {
            // Flatten per-sample SCO directories into [gene_id, fna_file] pairs,
            // then group by gene_id to get one tuple per gene across all samples.
            busco_gene_ch = BUSCO.out.single_copy_dir
                .flatMap { sample_id, scogs_dir ->
                    scogs_dir.listFiles()
                        .findAll { it.name.endsWith('.faa') }
                        .collect { f -> [f.baseName, f] }
                }
                .groupTuple()
                .filter { gene_id, files -> files.size() >= 2 }

            // Subsample genes for supermatrix if requested (default: 100)
            if (params.busco_subsample != 'all') {
                int n_genes = params.busco_subsample as Integer
                busco_gene_ch = busco_gene_ch
                    .toList()
                    .flatMap { all_genes ->
                        if (all_genes.size() <= n_genes) {
                            return all_genes
                        }
                        def shuffled = new ArrayList(all_genes)
                        Collections.shuffle(shuffled)
                        return shuffled.subList(0, n_genes)
                    }
            }

            // Align each gene across samples
            AlignBuscoGene(busco_gene_ch)

            // Build supermatrix (Python step, separate container from IQ-tree)
            concat_script_ch = Channel.value(file("${projectDir}/bin/concat_busco_alignments.py"))
            BuildSupermatrix(AlignBuscoGene.out.alignment.collect(), concat_script_ch)

            // Infer phylogeny with RAxML-NG (default) or IQ-tree
            if (params.busco_tree == 'iqtree') {
                IQtreePhylogeny(BuildSupermatrix.out.supermatrix, BuildSupermatrix.out.partitions)
                phylo_done = IQtreePhylogeny.out.results.collect()
            } else {
                RaxmlNgPhylogeny(BuildSupermatrix.out.supermatrix, BuildSupermatrix.out.partitions)
                phylo_done = RaxmlNgPhylogeny.out.results.collect()
            }
        }
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

        // [query_id, query_file] — used by both BLAST and Exonerate paths
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

        // ------------------------------------------------------------------
        // Paralog check (--query_paralog)
        // The paralogs of each query are BLASTed against every assembly; at the
        // hit locus of each sample, the query and paralog scores are compared.
        // ------------------------------------------------------------------
        if (params.query_paralog) {
            // One FASTA per query with all its paralogs ('<query>_paralogN'),
            // keyed by the safe form of the query's first header word
            paralog_fastas = channel.fromPath(params.query_paralog)
                .splitFasta(record: [id: true, text: true])
                .map { rec -> [rec.id.replaceAll(/_paralog[0-9]+$/, ''), rec.text] }
                .collectFile { key, text -> ["${safeName(key)}_paralogs.fasta", text] }
                .map { f -> [f.name.replaceAll(/_paralogs\.fasta$/, ''), f] }

            // [query_id, query_fasta, seq_type, paralog_fasta]
            paralog_queries = queries_with_type
                .map { query_id, query_fasta, seq_type ->
                    [safeName(query_id.trim().tokenize(' \t')[0]), query_id, query_fasta, seq_type]
                }
                .join(paralog_fastas)
                .map { _key, query_id, query_fasta, seq_type, paralog_fasta ->
                    [query_id, query_fasta, seq_type, paralog_fasta]
                }

            BLASTparalogs(
                all_scaffolds_ch
                    .combine(paralog_queries)
                    .map { sample_id, scaffolds_fasta, query_id, _query_fasta, seq_type, paralog_fasta ->
                        [sample_id, scaffolds_fasta, query_id, paralog_fasta, seq_type]
                    }
            )

            query_blast_by_query = BLASTsearch.out.results
                .map { query_id, _sample_id, results -> [query_id, results] }
                .groupTuple()
            paralog_blast_by_query = BLASTparalogs.out.results
                .map { query_id, _sample_id, results -> [query_id, results] }
                .groupTuple()

            ParalogCheck(
                paralog_queries
                    .map { query_id, query_fasta, _seq_type, paralog_fasta -> [query_id, query_fasta, paralog_fasta] }
                    .join(CombineBlastHits.out.summary)
                    .join(query_blast_by_query)
                    .join(paralog_blast_by_query)
            )
            paralog_done = ParalogCheck.out.table.map { _query_id, tsv -> tsv }.collect()
        }

        // Align BLAST hit sequences across samples
        AlignMAFFT(CombineBlastHits.out.combined_sequences.map { query_id, fasta -> [query_id, "blast", fasta] })

        // Build a phylogeny for each query from the BLAST hit alignment
        IQtreeBlastPhylogeny(
            AlignMAFFT.out.alignment.map { query_id, fasta -> [query_id, "blast", fasta] }
        )

        // ------------------------------------------------------------------
        // Exonerate path  (--search_tool exonerate)
        // Runs on BLAST-extracted regions (not the full assembly) for speed.
        // Performs protein-to-genome or cDNA-to-genome spliced alignment to
        // resolve exact exon boundaries. Outputs nucleotide CDS and protein.
        // ------------------------------------------------------------------
        if (params.search_tool == 'exonerate') {

            // Build Exonerate input from BLAST-extracted hit FASTAs (not full scaffolds):
            // [sample_id, hit_fasta, query_id, query_fasta, seq_type]
            exonerate_input = ExtractBlastHits.out.sequences
                .combine(queries_with_type, by: 0)
                .map { query_id, sample_id, hit_fasta, query_fasta, seq_type ->
                    [sample_id, hit_fasta, query_id, query_fasta, seq_type]
                }

            ExonerateSearch(exonerate_input)

            // Extract CDS + protein from the raw GFF3 using python3 (spades container)
            ExtractExonerateCds(ExonerateSearch.out.gff3)

            // Classify each gene model as full-length or truncated (nucleotide-level
            // deletion vs ORF mutation). The full assembly gives exact contig ends.
            integrity_input = ExonerateSearch.out.gff3
                .combine(queries_with_type, by: 0)
                .map { query_id, sample_id, hit_fasta, gff3, query_fasta, seq_type ->
                    [sample_id, query_id, hit_fasta, gff3, query_fasta, seq_type]
                }
                .combine(all_scaffolds_ch, by: 0)
                .map { sample_id, query_id, hit_fasta, gff3, query_fasta, seq_type, assembly ->
                    [query_id, sample_id, hit_fasta, gff3, query_fasta, seq_type, assembly]
                }
            AssessCdsIntegrity(integrity_input)

            // Group per-sample outputs by query_id
            exonerate_cds_by_query = ExtractExonerateCds.out.cds_sequences
                .map { query_id, sample_id, cds_fasta -> [query_id, cds_fasta] }
                .groupTuple()

            exonerate_prot_by_query = ExtractExonerateCds.out.protein_sequences
                .map { query_id, sample_id, prot_fasta -> [query_id, prot_fasta] }
                .groupTuple()

            exonerate_summary_by_query = ExtractExonerateCds.out.summary
                .groupTuple()

            // Join query file + all per-sample outputs for each query_id
            exonerate_combined_input = individual_queries
                .join(exonerate_cds_by_query)
                .join(exonerate_prot_by_query)
                .join(exonerate_summary_by_query)

            CombineExonerateHits(exonerate_combined_input)

            // Merge per-query exonerate summary TSVs into a single wide-format table
            all_exonerate_summaries = CombineExonerateHits.out.summary
                .map { query_id, summary_file -> summary_file }
                .collect()
            MergeExonerateSummaries(all_exonerate_summaries)

            // Align nucleotide CDS sequences across samples
            AlignExonerateCdsMAFFT(CombineExonerateHits.out.combined_cds.map { query_id, fasta -> [query_id, "exonerate", fasta] })

            // Align protein sequences across samples
            AlignProteinMAFFT(CombineExonerateHits.out.combined_protein.map { query_id, fasta -> [query_id, "exonerate", fasta] })

            // Build a phylogeny for each query from the CDS and protein alignments
            IQtreeCdsPhylogeny(
                AlignExonerateCdsMAFFT.out.alignment.map { query_id, fasta -> [query_id, "cds", fasta] }
            )
            IQtreeProteinPhylogeny(
                AlignProteinMAFFT.out.alignment.map { query_id, fasta -> [query_id, "protein", fasta] }
            )

            // Isoform summary: join CDS and protein alignments per query and assign
            // consecutive haplotype numbers for unique CDS / protein sequences.
            // The per-sample CDS integrity tables flag truncated isoforms; a query
            // without any (all AssessCdsIntegrity tasks failed) still gets a summary.
            integrity_by_query = AssessCdsIntegrity.out.summary
                .groupTuple()
                .join(AssessCdsIntegrity.out.events.groupTuple())
            isoform_input = AlignExonerateCdsMAFFT.out.alignment
                .join(AlignProteinMAFFT.out.alignment)
                .join(integrity_by_query, remainder: true)
                .filter { row -> row[1] != null }
                .map { row ->
                    def (query_id, cds_aln, prot_aln, summaries, events) = row + [null, null]
                    [query_id, cds_aln, prot_aln, summaries ?: [], events ?: []]
                }
            SummarizeIsoforms(isoform_input)

            // Isoform trees: one aligned row per CDS / protein isoform
            IQtreeCdsIsoformPhylogeny(
                SummarizeIsoforms.out.cds_isoforms_aligned.map { query_id, fasta -> [query_id, "cds_isoform", fasta] }
            )
            IQtreeProteinIsoformPhylogeny(
                SummarizeIsoforms.out.protein_isoforms_aligned.map { query_id, fasta -> [query_id, "protein_isoform", fasta] }
            )

            // Render each isoform tree as a PDF, labelled from the isoform FASTA
            // headers. Skipped when IQ-TREE wrote no tree (fewer than 3 isoforms).
            cds_isoform_trees = IQtreeCdsIsoformPhylogeny.out.results
                .filter { query_id, files -> (files instanceof List ? files : [files]).any { f -> f.name.endsWith('.treefile') } }
                .join(SummarizeIsoforms.out.cds_isoforms)
                .map { query_id, files, fasta -> [query_id, "cds_isoform", files, fasta] }
            protein_isoform_trees = IQtreeProteinIsoformPhylogeny.out.results
                .filter { query_id, files -> (files instanceof List ? files : [files]).any { f -> f.name.endsWith('.treefile') } }
                .join(SummarizeIsoforms.out.protein_isoforms)
                .map { query_id, files, fasta -> [query_id, "protein_isoform", files, fasta] }
            PlotIsoformTree(cds_isoform_trees.mix(protein_isoform_trees))

            // TEs on the scaffold carrying each Exonerate gene model (--repeat_library)
            if (params.repeat_library) {
                SummarizeHitTEs(CombineExonerateHits.out.summary.combine(repeat_out_all))
                hit_tes_done = SummarizeHitTEs.out.summary.map { query_id, tsv -> tsv }.collect()
            }

            exonerate_done = CombineExonerateHits.out.combined_cds
                .map { query_id, fasta -> fasta }

            // ------------------------------------------------------------------
            // Locus coverage (--locus_coverage)
            // ------------------------------------------------------------------
            if (params.locus_coverage) {
                blast_locus_input = queries_with_type
                    .map { query_id, query_fasta, seq_type ->
                        [query_id, query_fasta, seq_type, file(params.reference_genome)]
                    }
                BlastLocusOnReference(blast_locus_input)

                locus_cov_input = bam_ch_for_locus
                    .combine(BlastLocusOnReference.out.locus)
                    .map { sample_id, bam, bai, query_id, locus_tsv ->
                        [query_id, sample_id, bam, bai, locus_tsv]
                    }
                ComputeLocusCoverage(locus_cov_input)

                cov_by_query = ComputeLocusCoverage.out.coverage
                    .map { query_id, sample_id, tsv -> [query_id, tsv] }
                    .groupTuple()
                MergeLocusCoverage(cov_by_query)

                locus_cov_done = MergeLocusCoverage.out.merged
                    .map { query_id, f -> f }.collect()
            }

            all_done = assembly_qc_done
                .mix(MergeBlastSummaries.out)
                .mix(MergeExonerateSummaries.out.merged_summary.collect())
                .mix(exonerate_done)
                .mix(IQtreeBlastPhylogeny.out.results.collect())
                .mix(IQtreeCdsPhylogeny.out.results.collect())
                .mix(IQtreeProteinPhylogeny.out.results.collect())
                .mix(PlotIsoformTree.out.pdf.collect())
                .mix(hit_tes_done)
                .mix(paralog_done)
                .mix(fastp_summary_done)
                .mix(bwa_done)
                .mix(repeat_done)
                .mix(phylo_done)
                .mix(locus_cov_done)
                .collect()

        } else {
            // TEs on the scaffold carrying the best BLAST hit (--repeat_library)
            if (params.repeat_library) {
                SummarizeHitTEs(CombineBlastHits.out.summary.combine(repeat_out_all))
                hit_tes_done = SummarizeHitTEs.out.summary.map { query_id, tsv -> tsv }.collect()
            }

            // ------------------------------------------------------------------
            // Locus coverage (--locus_coverage)
            // ------------------------------------------------------------------
            if (params.locus_coverage) {
                blast_locus_input = queries_with_type
                    .map { query_id, query_fasta, seq_type ->
                        [query_id, query_fasta, seq_type, file(params.reference_genome)]
                    }
                BlastLocusOnReference(blast_locus_input)

                locus_cov_input = bam_ch_for_locus
                    .combine(BlastLocusOnReference.out.locus)
                    .map { sample_id, bam, bai, query_id, locus_tsv ->
                        [query_id, sample_id, bam, bai, locus_tsv]
                    }
                ComputeLocusCoverage(locus_cov_input)

                cov_by_query = ComputeLocusCoverage.out.coverage
                    .map { query_id, sample_id, tsv -> [query_id, tsv] }
                    .groupTuple()
                MergeLocusCoverage(cov_by_query)

                locus_cov_done = MergeLocusCoverage.out.merged
                    .map { query_id, f -> f }.collect()
            }

            all_done = assembly_qc_done
                .mix(MergeBlastSummaries.out)
                .mix(IQtreeBlastPhylogeny.out.results.collect())
                .mix(hit_tes_done)
                .mix(paralog_done)
                .mix(fastp_summary_done)
                .mix(bwa_done)
                .mix(repeat_done)
                .mix(phylo_done)
                .mix(locus_cov_done)
                .collect()
        }

    } else {
        all_done = assembly_qc_done
            .mix(fastp_summary_done)
            .mix(bwa_done)
            .mix(repeat_done)
            .mix(phylo_done)
            .mix(locus_cov_done)
            .collect()
    }

    // ---------------------
    // Generate pipeline execution statistics
    // ---------------------
    PipelineStatistics(all_done, workflow.launchDir.resolve(params.outdir).resolve("pipeline_trace.txt"))
}
