# genomedenovo_nf — Example dataset

## Source

All data are from the **E. coli Long-Term Evolution Experiment (LTEE)**,
BioProject [PRJNA295606](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA295606/).

> Tenaillon O, Barrick JE, Ribeck N, Deatherage DE, Blanchard JL, Dasgupta A,
> Wu GC, Wielgoss S, Schneider D, Blount ZD, Lenski RE (2016).
> Tempo and mode of genome evolution in a 50,000-generation experiment.
> *Nature* **536**, 165–170. https://doi.org/10.1038/nature18959

---

## Files

### Local FASTQ (pre-downloaded, ~19× coverage each)

Each file was streamed from ENA FTP and truncated to 300,000 read pairs
(150 bp PE) to limit file size for demonstration purposes.
Original full runs are available from ENA/SRA at the accessions below.

| File | SRR accession | Clone | Population | Generation |
|------|--------------|-------|------------|-----------|
| `fastq/SRR2589044_{1,2}.fastq.gz` | SRR2589044 | REL2181A | Ara−3 | 5,000 |
| `fastq/SRR2584863_{1,2}.fastq.gz` | SRR2584863 | REL7179B | Ara−3 | 15,000 |
| `fastq/SRR2584866_{1,2}.fastq.gz` | SRR2584866 | REL11365 | Ara−3 | 50,000 |

> **Note:** `fastq/` is excluded from git (see `.gitignore`).
> Re-download using the commands in the next section.

### Additional SRA accessions (downloaded at runtime)

Three further LTEE runs from different populations are listed in
`sra_accessions.txt` and will be fetched automatically when
`--SRA_index` is used.

| SRR accession | Clone | Population | Generation |
|--------------|-------|------------|-----------|
| SRR2591045 | REL7180A | Ara−4 | 15,000 |
| SRR2584867 | REL765A  | Ara−4 | 500 |
| SRR2591036 | REL4532A | Ara+3 | 10,000 |

### BLAST query sequence

| File | Accession | Gene | Description |
|------|-----------|------|-------------|
| `lacZ_WP_000177906.fasta` | WP_000177906.1 | *lacZ* | Beta-galactosidase protein from *E. coli* REL606 (NC_012967.1, complement 335787..338861). Retrieved via NCBI Entrez (`efetch -db protein -id WP_000177906.1 -format fasta`). |

---

## Running the example

Run from the repository root directory.

**Option 1 — local FASTQ files only (3 samples, assembly + BUSCO + lacZ BLAST):**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --blast_query example/lacZ_WP_000177906.fasta \
    --busco_lineage bacteria_odb10 \
    --outdir example/output
    --
```

**Option 2 — local FASTQ files + 3 additional SRA downloads (6 samples total):**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --SRA_index example/sra_accessions.txt \
    --blast_query example/lacZ_WP_000177906.fasta \
    --outdir example/output
```

**Option 3 - local FASTQ files + Exonerate search for lacZ (after BLAST):**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --blast_query example/lacZ_WP_000177906.fasta \
    --search_tool exonerate \
    --outdir example/output
```

**Option 4 - produce BUSCO gene phylogeny:**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --blast_query example/lacZ_WP_000177906.fasta \
    --busco_lineage bacteria_odb10 \
    --phylo_tree \
    --outdir example/output
```

**Option 5 - produce mapped read files (BAM) by providing a reference genome, assess read coverage of query hit locus:**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --blast_query example/lacZ_WP_000177906.fasta \
    --busco_lineage bacteria_odb10 \
    --reference_genome example/ecoli_K12_MG1655.fasta \
    --locus_coverage \
    --outdir example/output
```

**Option 6 - run repeat masking with a custom repeat library:**

```bash
nextflow run main.nf \
    --reads "example/fastq/SRR*_{1,2}.fastq.gz" \
    --blast_query example/lacZ_WP_000177906.fasta \
    --repeat_library example/ecoli_IS_elements.fasta \
    --outdir example/output
```