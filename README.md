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

A typical command for computing core genome variant calls from short-read data with Snippy:

    nextflow run np-core/np-variants --workflow core --fastq "isolates/*_R{1,2}.fq.gz" --reference ref.fasta

```

## Core Variants (Illumina)

In its simplest incarnation, the variant calling workflow uses `Snippy` and `Gubbins` to generate a non-recombinant core-genome alignment for the provided set of isolates (includes output of non-core per-isolate variant calls) from high-quality short-read sequence data:

```
nextflow run np-core/np-variants --workflow core --fastq "isolates/*_R{1,2}.fq.gz" --reference ref.fasta
```

Modules:

* `Fastp` - quality control of `--fastq` sequence reads 
* `Snippy` - reference alignment and core variant calls from `--fastq` and `--fasta`
* `Gubbins` - removal of recombinant sites from the reference variant alignment
* `PhyBeast` - support module to remove of all non-polymorphic sites after `Gubbins`

## Candidate Variants (ONT)

Candidate variants (usually single nucleotide polymorphisms, for example core genome variants called with `--workflow core`) can be used as input to anchor the neural network basecalling output from `Guppy` and call candidates with `Megalodon`. `NanoPath` then merges candidates and calls into a filtered alignment which can be used for 'hybrid' phylogenetic trees, combining both Illumina and ONT sequence data. While this is useful to reconsruct divergences within a known evolutionary background, for example when producing barcoded nanopore panels from an outbreak of a known lineage for which sufficient background data is available, please note that within-outbreak branches would not consider novel variation (particularly if divergence is deep, that is, if the outbreak has been persisting and accumulated novel variation for some time) and thus, within-oubtreak phylodynamic estimates based on candidate variants may not be accurate.

`Guppy` is used for basecalling and should be run using `GPU` resourcing through configuration files and profiles - if you require more information regarding the selection of resources for this workflow, please visit: [`np-core/configs`](https://github.com/np-core/configs)

```
nextflow run np-core/np-variants --config jcu -profile tesla --workflow candidate --fast5 fast5/ --candidates core.vcf
```

Modules used:

* `Megalodon` - basecalling and variant anchoring using `Guppy`
* `NanoPath` - merging and filtering of candidate and called variants

## *De novo* Variants (ONT)

*De novo* variant calling from nanopore sequence data and variant polishing

**under construction*

Variant polishers trained with `Nextflow` sub-workflow: `forest_training` and `snippy`

Basecall isolates with the same model (HAC methylation recommended) and call reference variants against the training reference genome with Snippy.

```
nextflow run np-core/np-signal --config jcu -profile tesla --caller guppy --fast5 fast5/ --candidates core.vcf
```

Create a directory holding one or more training directories (collections) with basecalled training isolates (`.fastq`) and corresponding reference variant (`.ref.vcf`), for example:

```
dar.fasta                  <-- param: --train_reference
train/                     <-- param: --dir_train
   saureus_mix/            <-- collection 1
      isolate1.fastq
      isolate1.ref.vcf
      isolate2.fastq
      isolate2.ref.vcf
   saureus_fnq/            <-- collection 2
      isolate3.fastq
      isolate3.ref.vcf
      isolate4.fastq
      isolate4.ref.vcf
```

Execute the training pipeline on the directory containing one or multiple training directories, for example:

```
nextflow run np-core/np-variants -profile docker --workflow forest_training --dir_train train/ --outdir dar_forest/ --train_reference dar.fasta
```


Modules used:

* `Clair` - haploid variant calling on human-trained models (default)
* `Medaka` - haploid variant calling on mixed-trained models (default)
