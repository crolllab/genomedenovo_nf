#!/usr/bin/env nextflow

params.outdir = "./nf_output/"
params.reads = ""
params.SRA_index = ""
params.NCBI_API_key = ""
params.blast_query = ""
params.scaffolds = ""
params.padding = 500  // Number of bp to extend BLAST hit sequences on both sides
params.SRR_sample_map = ""  // Optional: CSV file (no header) mapping original sample IDs to new names
params.search_tool = 'blast'  // Sequence search tool: 'blast' (default) or 'gth' (GenomeThreader)
params.skip_quast = false      // Set to true to skip QUAST assembly QC


include { denovo_query } from './workflows/denovo_query.nf'

log.info """

   ==========================================
   ||    GENOME DENOVO & QUERY WORKFLOW.   ||
   ==========================================

    Steps:
     1. SRA resolve & download (optional)
     2. De novo assembly (SPAdes)
     3. Assembly QC (QUAST)
     4. Sequence search: BLAST (always) + GenomeThreader CDS/protein (--search_tool gth)
     5. Hit extraction & per-query FASTA
     6. Multiple sequence alignment (MAFFT)
   
   Configuration:
     outdir         : ${params.outdir}
     local fastq    : ${params.reads}
     local scaffolds: ${params.scaffolds}
     blast query    : ${params.blast_query}
     SRA ids file   : ${params.SRA_index}
     NCBI API key   : ${params.NCBI_API_key}
     search tool    : ${params.search_tool}
     padding        : ${params.padding} bp
     sample map     : ${params.SRR_sample_map ?: 'none'}
    """


workflow {
   denovo_query()
}
