# Summary: BLAST Hit Padding Feature

**Added**: October 3, 2025  
**Parameter**: `--padding`  
**Default**: 0 (no padding)

## Quick Overview

The `--padding` parameter extends extracted BLAST hit sequences by adding flanking regions on both sides.

## Usage

```bash
# No padding (default) - exact BLAST hit only
nextflow run main.nf --blast_query genes.fa --scaffolds assemblies/

# Add 100 bp flanking on each side
nextflow run main.nf --blast_query genes.fa --scaffolds assemblies/ --padding 100

# Add 500 bp flanking (good for promoters)
nextflow run main.nf --blast_query genes.fa --scaffolds assemblies/ --padding 500
```

## What It Does

**Example**: BLAST hit at positions 1000-2000 with `--padding 100`

```
Without padding:     |====1000-2000====|  (1000 bp)
With padding 100:  |==900=========2100==|  (1200 bp)
```

- Extracts 900-2100 instead of 1000-2000
- Adds 100 bp upstream + hit + 100 bp downstream
- Automatically bounded to scaffold limits (won't go negative or beyond scaffold end)

## Output

### FASTA Headers Show Padded Coordinates
```
>sample1_wagA_hit_NODE_225:11901-15151
# Extracted region: 11901-15151 (with padding)
```

### TSV Files Show Original BLAST Coordinates
```tsv
Sample_ID  Query_ID  Subject_ID  ...  Start   End    ...
sample1    wagA      NODE_225    ...  12001   15051  ...
# BLAST alignment: 12001-15051 (original)
```

## Common Use Cases

| Use Case | Recommended Padding |
|----------|-------------------|
| Exact hit only | 0 (default) |
| Primer design | 100-200 bp |
| Promoter capture | 500-1000 bp |
| Gene cluster | 2000-5000 bp |
| Large context | 5000+ bp |

## Key Features

✅ Automatic boundary checking (won't exceed scaffold limits)  
✅ Works with forward and reverse strand hits  
✅ TSV preserves original BLAST coordinates  
✅ FASTA headers show extracted coordinates  
✅ No performance impact for reasonable padding values  

## Files Modified

- `main.nf` - Added parameter definition
- `modules/extract_blast_hits.nf` - Implemented padding logic
- `README.md` - Added parameter documentation
- `FEATURE_padding.md` - Complete technical documentation

## Example Output Comparison

### Default (no padding):
```bash
--padding 0  # Extract 1000-2000 (exact hit, 1000 bp)
```

### With padding:
```bash
--padding 100  # Extract 900-2100 (1200 bp total)
--padding 500  # Extract 500-2500 (2000 bp total)
```

Perfect for capturing regulatory regions, primers, or genomic context! 🎯
