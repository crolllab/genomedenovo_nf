// QUAST process for assembly quality assessment
process QUAST {
    tag "QUAST assessment: ${sample_id}"
    errorStrategy 'ignore'
    memory '4GB'
    cpus 2
    publishDir "${params.outdir}/quast_reports", mode: 'copy', pattern: "${sample_id}_quast"
    
    input:
    tuple val(sample_id), path(scaffolds_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_quast"), emit: report_dir
    path "${sample_id}_report.tsv", emit: metrics
    
    script:
    """
    quast.py ${scaffolds_fasta} \
        -o ${sample_id}_quast \
        --threads ${task.cpus} \
        --no-plots \
        --no-html \
        --no-icarus
    
    # Copy and rename the report.tsv with sample ID to avoid name collisions
    cp ${sample_id}_quast/report.tsv ${sample_id}_report.tsv
    
    # Add sample ID column to the TSV for easier merging
    sed -i "1s/Assembly/${sample_id}/" ${sample_id}_report.tsv
    """
}

// Process to collect and summarize all QUAST metrics
process SummarizeAssemblyMetrics {
    tag "Summarizing assembly metrics"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path quast_reports
    
    output:
    path "assembly_summary_metrics.tsv"
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os
    
    # Find all QUAST report TSV files (now named with sample IDs)
    tsv_files = glob.glob("*_report.tsv")
    
    if not tsv_files:
        # Create empty summary if no files found
        with open("assembly_summary_metrics.tsv", "w") as f:
            f.write("Sample_ID\\tNo_assembly_metrics_found\\n")
    else:
        # Read and combine all TSV files
        dataframes = []
        
        for file in tsv_files:
            try:
                # Extract sample name from filename
                sample_name = file.replace("_report.tsv", "")
                
                # Read the TSV (it's transposed format from QUAST)
                df = pd.read_csv(file, sep="\\t")
                
                # Transpose so metrics are columns
                df_transposed = df.set_index(df.columns[0]).T
                
                # Add sample ID as first column
                df_transposed.insert(0, "Sample_ID", sample_name)
                
                dataframes.append(df_transposed)
                
            except Exception as e:
                print(f"Error processing {file}: {e}")
        
        if dataframes:
            # Combine all dataframes
            combined_df = pd.concat(dataframes, ignore_index=True, sort=False)
            
            # Save combined metrics
            combined_df.to_csv("assembly_summary_metrics.tsv", sep="\\t", index=False)
        else:
            # Create empty summary if processing failed
            with open("assembly_summary_metrics.tsv", "w") as f:
                f.write("Sample_ID\\tError_processing_metrics\\n")
    """
}
