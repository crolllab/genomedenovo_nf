#!/usr/bin/env nextflow

params.outdir = "./nf_output/"
params.reference = ""
params.reads = ""
params.SRA_index = ""
params.ploidy = ""
params.NCBI_API_key = ""
params.keep_bam_gvcf = "false"  // or "true"
params.blast_query = ""
params.scaffolds = ""


include { denovo_query } from './workflows/denovo_query.nf'

log.info """

   ==========================================
   ||    GENOME DENOVO & QUERY WORKFLOW.   ||
   ==========================================

    Steps:
     1. SRA query and download
   
   Configuration:
     outdir         : ${params.outdir}
     reference      : ${params.reference}
     local fastq    : ${params.reads}
     local scaffolds: ${params.scaffolds}
     blast query    : ${params.blast_query}
     SRA ids file   : ${params.SRA_index}
     NCBI API key   : ${params.NCBI_API_key}
    """


workflow {
   denovo_query()
}
