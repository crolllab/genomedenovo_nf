// Separate processes for PE and SE downloads
process SRAdownloadPE {
    cpus 1
    memory '4GB'
    tag "Downloading $srr from NCBI SRA"
    
    input:
    val srr
    
    output:
    tuple val(srr), path("${srr}_1.fastq"), path("${srr}_2.fastq")
    
    script:
    """
    fasterq-dump $srr --split-files -O . 1>/dev/null
    
    # Verify we got paired files
    if [ ! -f "${srr}_2.fastq" ]; then
        echo "ERROR: Expected paired-end but only found single file for $srr" >&2
        exit 1
    fi
    """
}

process SRAdownloadSE {
    cpus 1
    memory '4GB'
    tag "Downloading $srr from NCBI SRA"
    
    input:
    val srr
    
    output:
    tuple val(srr), path("${srr}.fastq")
    
    script:
    """
    fasterq-dump $srr --split-files -O . 1>/dev/null
    
    # Verify we got single file only
    if [ -f "${srr}_2.fastq" ]; then
        echo "ERROR: Expected single-end but found paired files for $srr" >&2
        exit 1
    fi
    """
}