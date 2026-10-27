#!/usr/bin/env nextflow

params.outdir = "./nf_output/"
params.reads = ""
params.SRA_index = ""
params.NCBI_API_key = ""
params.blast_query = ""
params.scaffolds = ""
params.padding = 500  // Number of bp to extend BLAST hit sequences on both sides
params.SRR_sample_map = ""  // Optional: CSV file (no header) mapping original sample IDs to new names
params.search_tool = 'blast'  // Sequence search tool: 'blast' (default) or 'exonerate' (spliced alignment)
params.exonerate_refine = 'region' // Exonerate alignment refinement: 'region' (default), 'full' or 'none'
params.query_paralog = ""       // Optional: FASTA with paralogs of the queries, named '<query>_paralogN' (triggers the paralog check)
params.paralog_ambiguous_ratio = 0.95 // Paralog check: call 'ambiguous' when the second-best score reaches this fraction of the best
params.min_align_prop = 0.8     // Per-query trees: leave out sequences covering less than this proportion of the core alignment
params.skip_quast = false      // Set to true to skip QUAST assembly QC
params.busco_lineage = ""      // BUSCO lineage (e.g. dothideomycetes_odb10). If empty, BUSCO is skipped.
params.busco_db_path = ""      // Optional: path to pre-downloaded BUSCO databases (enables offline mode)
params.keep_full_busco = false // Optional: keep full BUSCO output directories (default: prune temporary BUSCO files)
params.reference_genome = ""   // Optional: reference genome FASTA for BWA-mem2 mapping (triggers bwa_mapping module)
params.repeat_library   = ""   // Optional: TE consensus library FASTA for RepeatMasker (triggers repeatmasker module)
params.phylo_tree       = false // Optional: build BUSCO-gene supermatrix phylogeny (requires --busco_lineage)
params.busco_tree       = 'raxml-ng' // Phylogenetic tool for BUSCO supermatrix: 'raxml-ng' (default) or 'iqtree'
params.busco_subsample  = 100   // Number of BUSCO genes to randomly subsample for the supermatrix (integer or 'all')
params.locus_coverage   = false // Optional: compute read depth at query BLAST locus on reference (requires --blast_query and --reference_genome)
params.locus_conflict_min_depth_ratio = 0.3 // Locus coverage: flag samples without the query in the assembly whose hit depth reaches this fraction of the genome depth
params.locus_conflict_min_breadth     = 80  // ... and whose hit breadth (% of bases with >= 1 read) reaches this value
params.slurm_queue      = ""    // Required when using -profile slurm: SLURM partition name (e.g. 'long', 'normal')


include { denovo_query } from './workflows/denovo_query.nf'

workflow {
   log.info """

   ==========================================
   ||    GENOME DENOVO & QUERY WORKFLOW.   ||
   ==========================================

    Steps:
     1. SRA resolve & download (optional)
     2. De novo assembly (SPAdes)
     3. Assembly QC (QUAST + BUSCO if --busco_lineage is set)
     4. Sequence search: BLAST (always) + Exonerate CDS/protein (--search_tool exonerate)
     5. Hit extraction & per-query FASTA
     6. Multiple sequence alignment (MAFFT) and per-query trees (IQ-TREE)
     7. Paralog check of the hits (--query_paralog)
   
   Configuration:
     outdir         : ${params.outdir}
     local fastq    : ${params.reads}
     local scaffolds: ${params.scaffolds}
     blast query    : ${params.blast_query}
     SRA ids file   : ${params.SRA_index}
     NCBI API key   : ${params.NCBI_API_key ? 'set' : 'not set'}
     search tool    : ${params.search_tool}
     exonerate refine: ${params.exonerate_refine}
     query paralogs : ${params.query_paralog ?: 'none (paralog check skipped)'}
     min align prop : ${params.min_align_prop}
     padding        : ${params.padding} bp
     sample map     : ${params.SRR_sample_map ?: 'none'}
     busco lineage  : ${params.busco_lineage ?: 'skipped'}
     keep full busco: ${params.keep_full_busco ? 'yes' : 'no (default compact outputs)'}
     reference      : ${params.reference_genome ?: 'none (bwa mapping skipped)'}
     repeat library : ${params.repeat_library ?: 'none (repeatmasker skipped)'}
     phylo tree     : ${params.phylo_tree ? "enabled (${params.busco_tree == 'iqtree' ? 'IQ-tree' : 'RAxML-NG'} BUSCO gene phylogeny, genes subsampled: ${params.busco_subsample})" : 'skipped'}
     SLURM queue    : ${params.slurm_queue ?: '(not set — not using SLURM)'}
    """

   if (workflow.profile.tokenize(',').contains('slurm') && !params.slurm_queue) {
       error """
    ERROR: --slurm_queue is required when using -profile slurm.
    Specify the SLURM partition name on your cluster, e.g.:
        --slurm_queue long
    The partition should allow a maximum walltime of at least 7 days for
    large datasets. Shorter limits may still work for smaller inputs, but
    jobs exceeding the partition limit will fail.
    """
   }

   denovo_query()
}
