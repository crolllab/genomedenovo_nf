![alt text](https://github.com/crolllab/genomedenovo_nf/blob/main/logo.png?raw=true)

Pipeline to assemble and query draft genomes
========================================

**genomedenovo_nf** is a [Nextflow](https://www.nextflow.io/) pipeline that assembles Illumina genomes de novo (or takes existing assemblies), checks assembly quality, and searches every genome for one or more query genes. For each query it extracts the hits, aligns them, builds trees and, optionally, recovers the exact coding sequences, classifies them as intact or truncated, and groups them into isoforms. It accepts paired-end reads from local files, NCBI/ENA SRA accessions, or pre-assembled scaffolds. Optional modules add a BUSCO species tree, a paralog check of the hits, read coverage at the query locus of a reference genome, and the transposable elements next to each hit.

Please report any issues or share feature requests on the [GitHub Issues](https://github.com/crolllab/genomedenovo_nf/issues) page.

## Installation and documentation

All documentation is here on the [genomedenovo_nf site](https://crolllab.github.io/genomedenovo_nf/).

See the [Getting started](https://crolllab.github.io/genomedenovo_nf/getting-started/) guide for installation instructions, example commands and input file formats. The [Configuration](https://crolllab.github.io/genomedenovo_nf/configuration/) page has a detailed description of all parameters and options. The [Output files](https://crolllab.github.io/genomedenovo_nf/output/) page describes the output files in detail, and the [Methods](https://crolllab.github.io/genomedenovo_nf/methods/) page what each step does. The [HPC & utilities](https://crolllab.github.io/genomedenovo_nf/resources/) page has recommendations for running on clusters. Release notes are in the [Changelog](https://crolllab.github.io/genomedenovo_nf/changelog/).
