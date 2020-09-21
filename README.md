# np-variants

![](https://img.shields.io/badge/lang-nextflow-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)
![](https://img.shields.io/badge/biorxiv-v0-f7fcb9.svg)

Reference alignment and variant calling for bacterial pathogens :orangutan:

```
nextflow np-core/np-variants --help true
```

## Usage

```
=========================================
   N P - V A R I A N T S  v${version}
=========================================

Usage:

A typical command for constructing the reference alignment and core-genome single nucleotide polymorphism calls:

    nextflow run np-core/np-variants --illumina fastq/ --tail "_R{1,2}.fastq.gz"

```

## Core Variants (Illumina)

In its simplest incarnation, the variant calling workflow uses `Snippy` and `Gubbins` to call non-recombinant core-genome variants for a set of isolates with Illumina PE reads.

```
nextflow run np-core/np-variants --fastq fastq/ --variants_sites true
```

Modules used:

* `Fastp` - quality control of `fastq` sequence reads 
* `Snippy` - reference alignment and core variant calls from `fastq` and `fasta`
* `Gubbins` - removal of recombinant sites from the reference variant alignment
* `Variant` - removal of all non-polymorphic sites after `Gubbins` (`--variant_sites true`)


## Megalodon Candidate Variants

`Megalodon` uses `Guppy` and should be run using `GPU` resourcing through configuration files and profiles - for more information see the `np-core/configs` repository.

```
nextflow run np-core/np-variants --config nextflow -profile gpu_docker --fast5 fast5/ --candidates core.vcf
```

Modules used:

* `Megalodon` which also uses `Guppy`

## ONT Variants

*De novo* variant calling from nanopore seqeunce data and assessment against reference `VCF`.
