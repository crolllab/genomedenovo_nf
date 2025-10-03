#!/usr/bin/env nextflow

params.outdir = "./nf_output/"
params.reads = ""
params.SRA_index = ""
params.NCBI_API_key = ""
params.blast_query = ""
params.scaffolds = ""
params.padding = 100  // Number of bp to extend BLAST hit sequences on both sides


include { denovo_query } from './workflows/denovo_query.nf'

log.info """

   ==========================================
   ||    GENOME DENOVO & QUERY WORKFLOW.   ||
   ==========================================

    Steps:
     1. SRA query and download
   
   Configuration:
     outdir         : ${params.outdir}
     local fastq    : ${params.reads}
     local scaffolds: ${params.scaffolds}
     blast query    : ${params.blast_query}
     SRA ids file   : ${params.SRA_index}
     NCBI API key   : ${params.NCBI_API_key}
     padding        : ${params.padding} bp
    """


workflow {
   denovo_query()
}
