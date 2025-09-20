#!/usr/bin/env nextflow

params.outdir = "./nf_output/"
params.reference = ""
params.reads = ""
params.SRA_index = ""
params.ploidy = ""
params.NCBI_API_key = ""
params.keep_bam_gvcf = "false"  // or "true"

include { variant_calling } from './workflows/variant_calling_wf'

log.info """

   ==========================================
   ||    GENOME DENOVO & QUERY WORKFLOW.   ||
   ==========================================

    Steps:
     1. SRA query and download
   
   Configuration:
     outdir         : ${params.outdir}
     reference      : ${params.reference}
     local files    : ${params.reads}
     SRA ids file   : ${params.SRA_index}
     NCBI API key   : ${params.NCBI_API_key}
    """


workflow {
   denovo_query()
}
