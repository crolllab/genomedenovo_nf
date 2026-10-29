# Configuration

All parameters are passed on the command line. Nextflow options take one dash (`-profile`), pipeline options two (`--reads`). Boolean flags are switched on by giving them (`--phylo_tree`) or with `true` / `false`. Optional steps are switched on by setting their option; everything else has a default.

---

## Startup validation

Parameters and input files are checked before any task is submitted. A missing input file, an invalid value (`--search_tool`, `--busco_tree`, `--busco_subsample`, `--min_align_prop`) or an option given without the options it depends on stops the run with an error message that names the problem. The table below lists the dependencies.

| Option | Requires |
|--------|----------|
| `--phylo_tree` | `--busco_lineage` |
| `--locus_coverage` | `--reference_genome`, `--blast_query` and read input (`--reads` or `--SRA_index`) |
| `--query_paralog` | `--blast_query` |
| `--search_tool exonerate` | Nucleotide queries without in-frame stop codons (a mature CDS) |
| `--busco_db_path` | The `--busco_lineage` folder inside `<busco_db_path>/lineages/` |
| `-profile slurm` | `--slurm_queue` |

The startup banner lists the resolved configuration, so you can check what the run will do before it starts. Sample names are checked for collisions once all inputs are known, before any assembly starts (see [Sample names](#sample-names)).

---

## File paths and glob patterns

`--reads` and `--scaffolds` take **glob patterns**; `--SRA_index`, `--blast_query`, `--query_paralog`, `--reference_genome`, `--repeat_library` and `--SRR_sample_map` take a single file.

### Quoting

**Always use single quotes around a glob pattern.**

| Form | Works? | What happens |
|------|--------|--------------|
| `--reads 'data/*_{1,2}.fastq.gz'` | ✅ **Use this** | The pattern reaches Nextflow intact and Nextflow expands it. |
| `--reads "data/*_{1,2}.fastq.gz"` | ⚠️ Usually | The shell does not glob-expand inside double quotes, so this normally behaves identically — but `$` and `` ` `` are still interpreted, so a path containing them breaks. |
| `--reads data/*_{1,2}.fastq.gz` | ❌ **Silently wrong** | The **shell** expands the pattern first. `--reads` receives only the *first* file and Nextflow discards the rest. |

No quoting is needed for a plain path with no wildcard in it, such as
`--reference_genome /data/genome.fasta`, though quoting does no harm.

### Reading from several locations

Give **one quoted value** and separate the patterns with a **semicolon**:

```bash
--reads 'runA/*_{1,2}.fastq.gz;runB/*_{1,2}.fastq.gz'
--scaffolds '/data/assemblies_2024;/data/assemblies_2025/*.fasta'
```

!!! warning "A comma is not a separator"
    A comma already means alternation *inside* a glob — the `{1,2}` in `*_{1,2}.fastq.gz` —
    so it cannot also separate patterns.

---

## Execution mode

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-profile` | `string` | none | Execution profile. `local` — local machine, 4 CPUs per task. `local_highCPU` — local machine, 16 CPUs per task. `slurm` — submit each task as a SLURM job (16 CPUs; requires `--slurm_queue`). Without a profile, tasks run locally with their own defaults (mostly 1 CPU). Tools whose CPU use was measured have their own CPU count in all profiles, see [Task resources](resources.md#task-resources). |
| `-resume` | flag | — | Reuse completed tasks from a previous run. Requires the `work-dir` to be intact. |
| `-work-dir` | `path` | `./work` | Directory for temporary and intermediate files, which can grow very large. Use a fast scratch filesystem for large datasets. Note the single dash. |
| `--outdir` | `path` | `./nf_output` | Directory for final output files. |
| `--slurm_queue` | `string` | — | **Required with `-profile slurm`.** Name of the SLURM partition to submit all jobs to. The partition should allow long jobs (up to 7 days for large assemblies). |
| `--max_cpus` | `integer` | `16` | Upper bound on the CPUs of any single task. |
| `--max_memory` | `memory` | `128.GB` | Upper bound on the memory of any single task, including memory that grows with each retry. |
| `--max_time` | `duration` | `7.d` | Upper bound on the walltime of any single task. |

!!! tip "Set the caps to your largest node"
    Requests above `--max_cpus`, `--max_memory` and `--max_time` are capped to these values instead of leaving a job pending forever. Set them to what one node (or your machine) offers.

---

## Read and assembly input

At least one of `--reads`, `--SRA_index` or `--scaffolds` is required; they can be combined.

### Local FASTQ files

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--reads` | `glob` | — | Local paired-end FASTQ files, e.g. `'/data/*_{1,2}.fq.gz'`. Single-quote the pattern; separate several locations with `;`. Reads are trimmed with fastp and assembled with SPAdes. |

### SRA / ENA accessions

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--SRA_index` | `path` | — | Plain-text file with one SRA/ENA accession per line (`SRR`, `SRX`, `SRP`, `PRJNA`, `ERR`, …). Every run of an accession is downloaded, from ENA when available and with NCBI `sra-tools` as fallback; paired-end and single-end runs are supported. [SRA Explorer](https://sra-explorer.info) helps to collect accessions. |
| `--NCBI_API_key` | `string` | — | NCBI API key for the SRA lookup, passed to `entrez-direct` as `NCBI_API_KEY`. Optional but recommended: it raises NCBI's request limit. Must be valid, or the lookup fails. Create one under *Account settings* of your [NCBI account](https://account.ncbi.nlm.nih.gov/). The startup banner only shows whether a key is set. |

### Existing assemblies

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--scaffolds` | `path` / `glob` | — | Existing assemblies (FASTA); trimming and assembly are skipped for them. Separate several sources with `;`. See [Scaffold input](#scaffold-input). |

#### Scaffold input

| `--scaffolds` value | Sample name |
|---|---|
| Folder with one subfolder per sample, each holding `scaffolds.fasta` (or `.fa`/`.fna`): `'/path/to/assemblies'` | Subfolder name |
| Recursive glob: `'/data/**/scaffolds.fasta'` | Parent folder name |
| File glob: `'/path/*.fasta'` | File name without extension (parent folder name if the file is called `scaffolds`, `assembly`, `genome` or `contigs`) |

### Sample names

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--SRR_sample_map` | `path` | — | CSV `original_id,new_name` (no header) to rename samples. Samples not listed keep their name. |

Sample names come from the input:

| Input | Sample name |
|---|---|
| `--reads` | Common prefix of the two read files, without trailing `_`, `.` or `-` (`S1_1.fq.gz` + `S1_2.fq.gz` → `S1`; note that `S1_R1`/`S1_R2` gives `S1_R`) |
| `--SRA_index` | SRR run accession |
| `--scaffolds` | Sample folder or file name (see [Scaffold input](#scaffold-input)) |

`--SRR_sample_map` renames them:

```
SRR1234567,Isolate_A
E1,Iran_2024_A3ghp1i
```

Inside the pipeline, and in all outputs, characters outside `[A-Za-z0-9._-]` are replaced by `_` (`Iran_2024_Kch(2)12-1g` → `Iran_2024_Kch_2_12-1g`), so that sample names are safe in file names, FASTA headers and trees. `<outdir>/sample_names.tsv` lists every input with its `--SRR_sample_map` name and final sample name.

!!! warning "Sample names must be unique"
    Final names must be unique (`isolate 1` and `isolate_1` clash). The pipeline checks this
    before assembly and stops with a list of the clashing inputs.

---

## Query search

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--blast_query` | `path` | — | FASTA with one or more query sequences, each searched separately. Nucleotide queries (≥ 95 % `ACGTUN`) use `blastn -task blastn`, protein queries `tblastn`. |
| `--search_tool` | `string` | `blast` | `blast`: BLAST only. `exonerate`: BLAST, then Exonerate spliced alignment for the exact CDS and protein, CDS integrity and isoforms. Nucleotide queries must then be a **mature CDS** (no introns); the run stops if a query has in-frame stop codons. |
| `--exonerate_refine` | `string` | `region` | Exonerate alignment refinement: `region` (re-align the hit region with full dynamic programming, so exons of divergent alleles are not lost), `full` (refine on the whole extracted window) or `none`. |
| `--padding` | `integer` | `500` | bp added on both sides of the BLAST hit when extracting it. The padded region is also the window searched by Exonerate and the locus used by `--locus_coverage`. With `--query_paralog`, the padding stops short of neighbouring paralog copies. |
| `--min_align_prop` | `float` | `0.8` | Per-query trees leave out sequences covering less than this proportion of the core alignment (the columns of the query reference, or of at least half of the sequences). The query reference is always kept. `0` keeps all sequences. See [Alignment filter](methods.md#alignment-filter). |

The per-query trees (BLAST hits, Exonerate CDS/protein, isoforms) do not need `--phylo_tree`: they are always built with `--blast_query`.

### Paralog check

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--query_paralog` | `path` | — | FASTA with known paralogs of the queries, switching on the [paralog check](methods.md#paralog-check). Name each paralog after the first word of its query header plus `_paralog1`, `_paralog2`, … (`>GeneA_paralog1` for the query `>GeneA …`); paralogs must be of the same type (nucleotide/protein) as their query. Queries without paralogs are not checked (a warning lists them). |
| `--paralog_ambiguous_ratio` | `float` | `0.95` | A query hit is called `ambiguous` when the second-best reference reaches this fraction of the query's bit score at the locus. Ambiguous hits are kept. |

!!! info "Paralog hits never enter the other outputs"
    With `--query_paralog`, each sample's query hit is extracted only where it belongs to the
    query. Samples whose query hits all lie on paralog loci are `paralog_only` in the
    summary tables and absent from the hit sequences, alignments and trees. See the
    [paralog tables](output.md#paralog-check-query_paralog).

---

## Assembly QC

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--skip_quast` | flag | off | Skip QUAST. |
| `--busco_lineage` | `string` | — | Run BUSCO with this OrthoDB lineage, e.g. `bacteria_odb10`, `fungi_odb10`, `ascomycota_odb10`, `dothideomycetes_odb10` ([all lineages](https://busco-data.ezlab.org/v5/data/lineages/)). Choose the most specific lineage for your organism. |
| `--busco_db_path` | `path` | — | Folder with pre-downloaded lineages (`busco --download <lineage> --download_path <folder>`); BUSCO then runs offline. Without it, the lineage is downloaded once into `<outdir>/busco_db_cache/`. |
| `--keep_full_busco` | flag | off | Keep the complete BUSCO output folders. By default only the short summaries and the single-copy protein sequences are kept. |

---

## BUSCO species tree

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--phylo_tree` | flag | off | Build a species tree from the BUSCO single-copy genes. Requires `--busco_lineage`. |
| `--busco_tree` | `string` | `raxml-ng` | `raxml-ng` or `iqtree`. |
| `--busco_subsample` | `integer` / `all` | `100` | Number of BUSCO genes randomly drawn for the supermatrix, or `all`. |

---

## Reference mapping and locus coverage

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--reference_genome` | `path` | — | Map the trimmed reads to this reference with BWA-mem2. Needs read input. |
| `--locus_coverage` | flag | off | Compare read depth at the query locus of the reference with the genome-wide depth. Needs `--reference_genome`, `--blast_query` and read input. |
| `--locus_conflict_min_depth_ratio` | `float` | `0.3` | A sample whose assembly lacks the query (no hit, or only paralog hits) is flagged when the read depth over the reference hit reaches this fraction of its genome depth… |
| `--locus_conflict_min_breadth` | `float` | `80` | …and at least this percentage of the hit's bases is covered by reads. See `coverage_conflict` in the [locus coverage table](output.md#11_locus_coverage). |

---

## Repeats

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--repeat_library` | `path` | — | TE consensus library (FASTA, e.g. from RepeatModeler or [Dfam](https://www.dfam.org)). All assemblies are annotated with RepeatMasker; with `--blast_query`, the TEs near each query hit are reported. |
