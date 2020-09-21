# np-variants

![](https://img.shields.io/badge/lang-nextflow-41ab5d.svg)
![](https://img.shields.io/badge/version-0.1.0-addd8e.svg)
![](https://img.shields.io/badge/biorxiv-v0-f7fcb9.svg)

Reference alignment and haploid variant calling for bacterial phylogenetics and -dynamics :orangutan:

```
nextflow np-core/np-variants --help true
```

- [Usage](#usage)
- Workflows
   - [Core Variants (Illumina)](#core-variants-illumina)
   - [Candiate Variants (ONT)](#candidate-variants-ont)
   - [*De novo* Variants (ONT)](#de-novo-variants-ont)

## Usage

```
=========================================
   N P - V A R I A N T S  v${version}
=========================================

Usage:

A typical command for constructing the reference alignment and core-genome single nucleotide polymorphism calls:

    nextflow run np-core/np-variants --workflow core --fastq isolates/ --tail "_R{1,2}.fastq.gz"

```

## Core Variants (Illumina)

In its simplest incarnation, the variant calling workflow uses `Snippy` and `Gubbins` to generate a non-recombinant core-genome alignment for the provided set of isolates (includes output of non-core per-isolate variant calls) from high-quality short-read sequence data:

```
nextflow run np-core/np-variants --workflow core --fastq isolates/ --fasta isolate_asemblies/ --variant_sites true
```

Modules used:

* `Fastp` - quality control of `--fastq` sequence reads 
* `Snippy` - reference alignment and core variant calls from `--fastq` and `--fasta`
* `Gubbins` - removal of recombinant sites from the reference variant alignment
* `PhyBeast` - support module to remove of all non-polymorphic sites after `Gubbins`

## Candidate Variants (ONT)

Candidate variants (usually single nucleotide polymorphisms, for example core genome variants called with `--workflow core`) can be used as input to anchor the neural network basecalling output from `Guppy` and call candidates with `Megalodon`. `NanoPath` then merges candidates and calls into a filtered alignment which can be used for 'hybrid' phylogenetic trees, combining both Illumina and ONT sequence data. While this is useful to reconsruct divergences within a known evolutionary background, for example when producing barcoded nanopore panels from an outbreak of a known lineage for which sufficient background data is available, please note that within-outbreak branches would not consider novel variation (particularly if divergence is deep, that is, if the outbreak has been persisting and accumulated novel variation for some time) and thus, within-oubtreak phylodynamic estimates based on candidate variants may not be accurate.

`Guppy` is used for basecalling and should be run using `GPU` resourcing through configuration files and profiles - if you require more information regarding the selection of resources for this workflow, please visit the [`np-core/configs`](https://github.com/np-core/configs) repository. In this example the local `JCU` configuration for our `Tesla` GPU server is used:

```
nextflow run np-core/np-variants --config jcu -profile tesla --workflow candidate --fast5 fast5/ --candidates core.vcf
```

Modules used:

* `Megalodon` - basecalling and variant anchoring using `Guppy`
* `NanoPath` - merging and filtering of candidate and called variants

## *De novo* Variants (ONT)

*De novo* variant calling from nanopore seqeunce data and assessment against reference `VCF`.
