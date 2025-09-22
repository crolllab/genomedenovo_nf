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
    import csv
    import glob
    import os
    
    # Find all QUAST report TSV files (now named with sample IDs)
    tsv_files = glob.glob("*_report.tsv")
    
    if not tsv_files:
        # Create empty summary if no files found
        with open("assembly_summary_metrics.tsv", "w") as f:
            f.write("Sample_ID\\tNo_assembly_metrics_found\\n")
    else:
        all_data = []
        all_metrics = set()
        
        # Process each TSV file
        for file in tsv_files:
            try:
                # Extract sample name from filename
                sample_name = file.replace("_report.tsv", "")
                
                # Read the TSV file
                with open(file, 'r') as f:
                    reader = csv.reader(f, delimiter='\\t')
                    rows = list(reader)
                
                if len(rows) < 2:
                    continue
                
                # QUAST format: first column is metric names, second column is values
                sample_data = {"Sample_ID": sample_name}
                
                for row in rows[1:]:  # Skip header
                    if len(row) >= 2:
                        metric_name = row[0].strip()
                        metric_value = row[1].strip()
                        sample_data[metric_name] = metric_value
                        all_metrics.add(metric_name)
                
                all_data.append(sample_data)
                
            except Exception as e:
                print(f"Error processing {file}: {e}")
        
        if all_data:
            # Sort metrics for consistent column order
            sorted_metrics = sorted(all_metrics)
            header = ["Sample_ID"] + sorted_metrics
            
            # Write combined TSV
            with open("assembly_summary_metrics.tsv", "w", newline='') as f:
                writer = csv.writer(f, delimiter='\\t')
                
                # Write header
                writer.writerow(header)
                
                # Write data rows
                for sample_data in all_data:
                    row = [sample_data.get(col, "") for col in header]
                    writer.writerow(row)
        else:
            # Create empty summary if processing failed
            with open("assembly_summary_metrics.tsv", "w") as f:
                f.write("Sample_ID\\tError_processing_metrics\\n")
    """
}
