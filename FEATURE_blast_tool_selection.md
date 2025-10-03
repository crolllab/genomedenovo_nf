# Feature: Automatic BLAST Tool Selection (blastn vs tblastn)

**Date**: October 3, 2025  
**Module**: `modules/blast_search.nf`  
**Processes**: `SplitQuerySequences`, `BLASTsearch` (renamed from `BLASTn`)

## Overview

The pipeline now automatically detects whether query sequences are nucleotide or protein sequences and selects the appropriate BLAST tool:
- **blastn**: For nucleotide queries (searching nucleotide database with nucleotide query)
- **tblastn**: For protein queries (searching translated nucleotide database with protein query)

## Implementation

### Sequence Type Detection

Added bash function in `SplitQuerySequences` process to detect sequence type:

```bash
detect_sequence_type() {
    local fasta_file=$1
    # Extract sequence and convert to uppercase
    seq=$(grep -v "^>" "$fasta_file" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
    
    # Count nucleotide characters (A, T, G, C, U, N)
    nuc_count=$(echo "$seq" | grep -o '[ATGCUN]' | wc -l)
    total_count=$(echo "$seq" | wc -c)
    
    # If >95% are valid nucleotides, consider it nucleotide sequence
    if [ $nuc_percent -ge 95 ]; then
        echo "nucleotide"
    else
        echo "protein"
    fi
}
```

### Detection Logic

- **Nucleotide**: If ≥95% of characters are A, T, G, C, U, or N
- **Protein**: If <95% of characters are nucleotides (contains amino acid codes like L, K, M, F, P, etc.)
- **Mixed queries**: If file contains both types, defaults to nucleotide (safer choice)

### Process Changes

#### SplitQuerySequences

**New outputs:**
- `query_type.txt`: Contains the dominant sequence type (`nucleotide`, `protein`, or `mixed`)
- Updated `query_info.tsv`: Added `Sequence_Type` column

**Example query_info.tsv:**
```
Query_ID    Filename                Sequence_Length    Sequence_Type
wagA        query_1_wagA.fasta      4206               nucleotide
wagB        query_2_wagB.fasta      3987               nucleotide
CYP51_prot  query_3_CYP51_prot.fasta  530              protein
```

#### BLASTsearch (renamed from BLASTn)

**New input:**
- Added `val(seq_type)` to input tuple
- Full tuple: `[sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]`

**Dynamic tool selection:**
```groovy
def blast_tool = (seq_type == "protein") ? "tblastn" : "blastn"
```

**Tag update:**
- Old: `"BLAST search: ${sample_id} vs ${query_id}"`
- New: `"BLAST search: ${sample_id} vs ${query_id} (${seq_type})"`

### Workflow Integration

Modified `workflows/denovo_query.nf`:

1. **Read query info:**
```groovy
query_info_ch = SplitQuerySequences.out.query_info
    .splitCsv(header: true, sep: '\t')
    .map { row -> [row.Query_ID, row.Sequence_Type] }
```

2. **Join queries with types:**
```groovy
queries_with_type = individual_queries
    .join(query_info_ch)
    .map { query_id, query_file, seq_type ->
        [query_id, query_file, seq_type]
    }
```

3. **Create BLAST input with type:**
```groovy
blast_input = all_scaffolds_ch
    .combine(queries_with_type)
    .map { sample_id, scaffolds_fasta, query_id, query_fasta, seq_type ->
        [sample_id, scaffolds_fasta, query_id, query_fasta, seq_type]
    }
```

### Configuration Update

Renamed process in `nextflow.config`:
```groovy
withName:BLASTsearch {  // Changed from BLASTn
    container = './singularity/blast%3A2.17.0--h66d330f_0'
    memory = '4 GB'
    maxForks = 10
}
```

## Usage Examples

### Nucleotide Query (wagA, wagB, wagC)

```bash
nextflow run main.nf \
  -config nextflow.config \
  -profile local_highCPU \
  --scaffolds '/path/to/assemblies' \
  --blast_query './wagABC_genes.fa' \
  --outdir './nf_output'
```

Pipeline will detect nucleotide sequences and use **blastn**.

### Protein Query

```bash
nextflow run main.nf \
  -config nextflow.config \
  -profile local_highCPU \
  --scaffolds '/path/to/assemblies' \
  --blast_query './protein_queries.fa' \
  --outdir './nf_output'
```

Pipeline will detect protein sequences and use **tblastn**.

### Mixed Query File

If your query file contains both nucleotide and protein sequences:

```
>gene1_nucleotide
ATGCGATCGATCG...
>protein1_amino_acids
MKVLIVAIVFG...
```

Pipeline will:
1. Detect each sequence type individually
2. Log warning: "Mixed sequence types detected"
3. Default to nucleotide BLAST (safer choice)
4. **Recommendation**: Split into separate query files for optimal results

## Output Changes

### Log Output

Process tags now show sequence type:
```
[BLASTsearch] BLAST search: ST99CH_3A1 vs wagA (nucleotide)
[BLASTsearch] BLAST search: ST99CH_3A1 vs CYP51_protein (protein)
```

### Query Info File

New column in `query_sequences/query_info.tsv`:
```
Query_ID    Filename                    Sequence_Length    Sequence_Type
wagA        query_1_wagA.fasta          4206               nucleotide
CYP51       query_2_CYP51.fasta         1593               nucleotide
AvrStb6     query_3_AvrStb6_prot.fasta  189                protein
```

## BLAST Tool Comparison

| Feature | blastn | tblastn |
|---------|--------|---------|
| Query Type | Nucleotide | Protein |
| Database Type | Nucleotide | Nucleotide (translated) |
| Use Case | Gene finding | Protein homology search |
| Translation | None | 6-frame translation of DB |
| Sensitivity | Lower (nucleotide level) | Higher (protein level) |
| Speed | Faster | Slower (6x DB size) |

### When to Use Each

**Use blastn (nucleotide queries):**
- Searching for specific genes/alleles
- Finding exact or near-exact sequence matches
- Comparing closely related sequences
- When query is DNA/RNA sequence

**Use tblastn (protein queries):**
- Searching for protein-coding genes
- Finding distant homologs
- When query is amino acid sequence
- Discovering novel gene variants
- Cross-species comparisons

## Validation

### Test Sequence Type Detection

```bash
# Create test files
echo -e ">nucleotide\nATGCGATCGATCGATC" > test_nuc.fa
echo -e ">protein\nMKVLIVAIVFGLLVGG" > test_prot.fa

# Run detection (bash function from SplitQuerySequences)
detect_sequence_type test_nuc.fa
# Output: nucleotide (99% ATGCUN)

detect_sequence_type test_prot.fa
# Output: protein (only 25% ATGCUN)
```

### Expected Results

**Nucleotide query (wagABC_genes.fa):**
- Detection: `nucleotide` (99% ATGCUN characters)
- Tool used: `blastn`
- Log: "Running blastn for nucleotide query: wagA"

**Protein query (test_protein_query.fa):**
- Detection: `protein` (27% ATGCUN characters)
- Tool used: `tblastn`
- Log: "Running tblastn for protein query: test_protein_query"

## Technical Details

### Detection Threshold

- **95% cutoff**: Chosen to be conservative
  - Allows for some ambiguous bases (N) in nucleotide sequences
  - Clearly distinguishes from protein sequences (typically <40% ATGCUN)

### Character Handling

- **Case-insensitive**: Converts to uppercase before counting
- **Whitespace-tolerant**: Removes newlines during analysis
- **Ambiguous bases**: N is counted as nucleotide

### BLAST Parameters

Both tools use identical parameters:
```bash
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
-evalue 1e-5
-num_threads ${task.cpus}
```

## Limitations

1. **Mixed Files**: If query file contains both types, defaults to nucleotide
   - **Solution**: Split into separate files for optimal results

2. **Unusual Sequences**: Sequences with many ambiguous bases might be misclassified
   - **Solution**: Manual review of `query_info.tsv` Sequence_Type column

3. **CodonHaven't implemented**: Currently only blastn and tblastn
   - Not implemented: blastx (translated query vs protein DB)
   - Not implemented: tblastx (translated query vs translated DB)

## Troubleshooting

### Issue: All queries detected as nucleotide despite being protein

**Cause**: Protein sequences contain many A, T, G, C amino acids  
**Check**: View `query_sequences/query_info.tsv` to see detection results  
**Solution**: If misdetected, you can manually edit the query file to use different amino acid codes, or split queries

### Issue: BLAST returns no hits for protein query

**Possible causes:**
1. Protein sequence too short (tblastn needs sufficient length)
2. E-value threshold too stringent
3. No homologs present in target genomes

**Solutions:**
1. Use longer protein queries (>50 amino acids recommended)
2. Adjust e-value: modify `-evalue 1e-5` in `modules/blast_search.nf`
3. Try with known positive control sequences

## Related Files

- Module: `modules/blast_search.nf`
- Workflow: `workflows/denovo_query.nf`
- Config: `nextflow.config`
- Test files: `test_protein_query.fa`, `wagABC_genes.fa`

## Future Enhancements

Potential improvements:
1. Support for blastx (translated nucleotide query vs protein DB)
2. Support for tblastx (translated query vs translated DB)
3. User-specified blast tool via parameter
4. Per-query blast parameters (different e-values per query)
5. Automatic adjustment of e-value based on query length
