# Feature: Padding for BLAST Hit Extraction

**Date**: October 3, 2025  
**Module**: `modules/extract_blast_hits.nf`  
**Process**: `ExtractBlastHits`  
**Parameter**: `--padding`

## Overview

Added a `--padding` parameter that extends extracted BLAST hit sequences by a specified number of base pairs on both upstream and downstream sides. This is useful for capturing flanking regions, regulatory elements, or providing additional context around BLAST hits.

## Implementation

### Parameter Definition

**Location**: `main.nf`
```groovy
params.padding = 0  // Number of bp to extend BLAST hit sequences on both sides
```

**Default**: 0 (extract exact BLAST hit region only)  
**Type**: Integer (number of base pairs)

### Usage

```bash
# Extract exact BLAST hit regions (default)
nextflow run main.nf --blast_query query.fa --scaffolds assemblies/

# Extract hits plus 100 bp flanking on each side
nextflow run main.nf --blast_query query.fa --scaffolds assemblies/ --padding 100

# Extract hits plus 500 bp flanking (useful for promoters)
nextflow run main.nf --blast_query query.fa --scaffolds assemblies/ --padding 500

# Extract hits plus 2000 bp flanking (large context)
nextflow run main.nf --blast_query query.fa --scaffolds assemblies/ --padding 2000
```

## Technical Details

### Coordinate Calculation

For each BLAST hit with coordinates `[start, end]`:

1. **Calculate padded coordinates**:
   ```bash
   padded_start = start - padding
   padded_end = end + padding
   ```

2. **Apply boundary constraints**:
   ```bash
   if padded_start < 1:
       padded_start = 1
   
   if padded_end > scaffold_length:
       padded_end = scaffold_length
   ```

3. **Extract sequence**:
   - Use `padded_start` and `padded_end` for extraction
   - Coordinates automatically stay within scaffold boundaries

### Example Calculation

**BLAST hit**: Coordinates 1000-2000 on a 10000 bp scaffold  
**Padding**: 100 bp

```
Original:      |----1000======2000----|  (1000 bp hit)
Padded:     |--900============2100--|  (1200 bp total: hit + 2×100)
```

**BLAST hit**: Coordinates 50-500 on a 1000 bp scaffold  
**Padding**: 100 bp

```
Scaffold:  |1===================1000|
Original:     |--50====500----|
Padded:      |1======600---|  (Start bounded at 1, not -50)
```

### File Output

#### FASTA Headers

FASTA sequence headers show **padded coordinates**:

```
>sample1_wagA_hit_NODE_225:11901-15151
# Original BLAST hit: 12001-15051
# With padding=100:    11901-15151 (added 100 bp on each side)
```

#### TSV Files

TSV hit info files contain **original BLAST coordinates** (not padded):

```tsv
Sample_ID    Query_ID  Subject_ID    Percent_Identity  Alignment_Length  Start  End    E_value  Bit_Score
sample1      wagA      NODE_225...   98.170           2951              12001  15051  0.0      5134
```

**Rationale**: 
- FASTA headers show what sequence was extracted (padded)
- TSV shows actual BLAST alignment (original)
- This allows users to know both the extraction region and alignment region

### Boundary Handling

#### Case 1: Padding Fits Within Scaffold

```
Scaffold: |1==================10000|
Hit:              |--5000===6000--|
Padding=500:    |4500===========6500|
Result: Extract 4500-6500 (2500 bp total)
```

#### Case 2: Padding Exceeds Start Boundary

```
Scaffold: |1==================10000|
Hit:      |50====500--|
Padding=100: |-50====600|  → Bounded to |1====600|
Result: Extract 1-600 (600 bp total, not 650)
```

#### Case 3: Padding Exceeds End Boundary

```
Scaffold: |1==================10000|
Hit:                      |9500===9900--|
Padding=200:            |9300=======10100| → Bounded to |9300====10000|
Result: Extract 9300-10000 (700 bp total, not 800)
```

#### Case 4: Hit Near Both Boundaries

```
Scaffold: |1==================10000|
Hit:      |10======100--|.........|9900==9990--|
Padding=50: |-40===150|           |9850==10040| 
Bounded:    |1====150|            |9850===10000|
```

### Reverse Strand Handling

Padding is applied to genomic coordinates **before** reverse complementing:

```
Forward strand hit: 1000-2000, padding=100
  → Extract 900-2100
  → Return sequence as-is

Reverse strand hit: 2000-1000 (BLAST reports reversed), padding=100
  → Convert to forward: 1000-2000
  → Apply padding: 900-2100
  → Extract sequence
  → Reverse complement
  → Return RC sequence
```

## Use Cases

### 1. Regulatory Region Analysis

**Goal**: Capture promoter regions upstream of genes

```bash
# Extract gene hits plus 500 bp upstream promoter
nextflow run main.nf \
  --blast_query genes.fa \
  --scaffolds assemblies/ \
  --padding 500
```

**Example**: Gene at positions 5000-6000
- Without padding: Extract 5000-6000 (gene only)
- With --padding 500: Extract 4500-6500 (promoter + gene + downstream)

### 2. Primer Design

**Goal**: Get flanking sequence for designing primers

```bash
# Extract hit plus 100 bp flanking for primer design
nextflow run main.nf \
  --blast_query target_regions.fa \
  --scaffolds assemblies/ \
  --padding 100
```

**Benefit**: 100 bp flanking provides sufficient sequence for specific primer design

### 3. Structural Variant Detection

**Goal**: Look for insertions/deletions in flanking regions

```bash
# Extract large flanking regions
nextflow run main.nf \
  --blast_query conserved_genes.fa \
  --scaffolds assemblies/ \
  --padding 2000
```

**Benefit**: Large padding allows detection of structural variants near hits

### 4. Gene Cluster Analysis

**Goal**: Capture neighboring genes

```bash
# Extract hit plus 5kb on each side for gene cluster
nextflow run main.nf \
  --blast_query cluster_marker.fa \
  --scaffolds assemblies/ \
  --padding 5000
```

### 5. Exact Hit Extraction (Default)

**Goal**: Extract only the aligned region

```bash
# No padding - exact BLAST hit only
nextflow run main.nf \
  --blast_query query.fa \
  --scaffolds assemblies/ \
  --padding 0  # or omit (default)
```

## Output Comparison

### Without Padding (--padding 0, default)

**FASTA**:
```
>sample1_wagA_hit_NODE_225:12001-15051
ATGCGATCGATC...
# Length: 3051 bp (exact BLAST hit)
```

**TSV**:
```
sample1  wagA  NODE_225  98.17  2951  12001  15051  0.0  5134
```

### With Padding (--padding 100)

**FASTA**:
```
>sample1_wagA_hit_NODE_225:11901-15151
ATGCGATCGATC...
# Length: 3251 bp (hit + 100 bp flanking on each side)
```

**TSV**:
```
sample1  wagA  NODE_225  98.17  2951  12001  15051  0.0  5134
# Note: TSV coordinates unchanged - show actual BLAST alignment
```

## Performance Considerations

### Memory Usage

Padding increases extracted sequence length:
- Padding = 0: ~3 KB per hit (typical 3000 bp gene)
- Padding = 100: ~3.2 KB per hit
- Padding = 500: ~4 KB per hit
- Padding = 2000: ~7 KB per hit

**Impact**: Minimal for typical use cases (<1000 samples)

### Execution Time

Padding does not significantly affect runtime:
- Coordinate calculation: Negligible (~0.001s per hit)
- Sequence extraction: Scales linearly with length
- For 1000 samples with padding=500: <2 seconds additional time

### Disk Space

Larger extracted sequences increase output file sizes:

**Example with 1000 samples, 3 queries, padding=500**:
- Without padding: ~9 MB combined FASTA
- With padding: ~12 MB combined FASTA (+33%)

## Edge Cases

### 1. Very Short Scaffolds

**Scenario**: Hit on 200 bp scaffold, padding=500

```
Scaffold: |1======200|
Hit:         |50===150|
Padding=500: |-450====650| → Bounded to |1====200|
Result: Extract entire scaffold (200 bp)
```

### 2. Hit Spans Entire Scaffold

**Scenario**: Hit covers whole scaffold

```
Scaffold: |1======1000|
Hit:      |1======1000|
Padding=100: |-100====1100| → Bounded to |1====1000|
Result: No change, already full scaffold
```

### 3. Multiple Hits on Same Scaffold

**Scenario**: Two hits on same scaffold, may overlap after padding

```
Scaffold: |1======================10000|
Hit A:    |1000==2000|
Hit B:              |2500==3500|
Padding=500:
  Hit A padded: |500========2500|
  Hit B padded:         |2000========4000|
Result: Two separate sequences extracted (overlaps OK)
```

Note: Overlapping extractions are allowed - each hit is processed independently.

### 4. Zero-Length or Invalid Hits

BLAST hits with invalid coordinates are automatically skipped (existing error handling).

## Validation

### Test Cases

```bash
# Test 1: No padding (default)
nextflow run main.nf --blast_query test.fa --scaffolds test/ --padding 0

# Test 2: Small padding
nextflow run main.nf --blast_query test.fa --scaffolds test/ --padding 50

# Test 3: Medium padding
nextflow run main.nf --blast_query test.fa --scaffolds test/ --padding 200

# Test 4: Large padding
nextflow run main.nf --blast_query test.fa --scaffolds test/ --padding 1000

# Test 5: Extreme padding (should be bounded)
nextflow run main.nf --blast_query test.fa --scaffolds test/ --padding 100000
```

### Verification

1. Check FASTA header coordinates match expected padding
2. Verify TSV coordinates are unchanged (original BLAST hits)
3. Ensure padded coordinates don't exceed scaffold boundaries
4. Confirm sequence length = (end - start + 1)

## Troubleshooting

### Issue: Extracted sequence shorter than expected

**Cause**: Padding extended beyond scaffold boundary  
**Check**: Compare FASTA header coordinates with scaffold length  
**Solution**: Normal behavior - coordinates are bounded to scaffold limits

### Issue: TSV coordinates don't match FASTA header

**Expected**: TSV shows original BLAST hit, FASTA shows padded extraction region  
**Solution**: This is intentional - both pieces of information are useful

### Issue: Overlapping sequences from same scaffold

**Cause**: Multiple hits with padding overlap  
**Solution**: Normal behavior - each hit extracted independently

## Future Enhancements

Potential improvements:
1. **Asymmetric padding**: Different upstream/downstream values
   - `--padding_upstream 500 --padding_downstream 100`
2. **Query-specific padding**: Different padding per query
   - Specified in query metadata file
3. **Merge overlapping hits**: Combine overlaps into single sequence
4. **Smart padding**: Extend to gene boundaries automatically

## Related Files

- **Modified**: `main.nf` (parameter definition)
- **Modified**: `modules/extract_blast_hits.nf` (extraction logic)
- **Updated docs**: `README.md`

## References

- samtools faidx documentation: http://www.htslib.org/doc/samtools-faidx.html
- BLAST output format: https://www.ncbi.nlm.nih.gov/books/NBK279684/
