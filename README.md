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
   - [Random Forest SNP polishers (ONT)](#de-novo-variants-ont)

## Usage

```
=========================================
        N P - V A R I A N T S  
=========================================

Usage:

A typical command to construct a bacterial core-genome SNP alignment from Illumina PE reads:

  nextflow run np-core/np-variants --workflow snippy --snippy_dir "*_R{1,2}.fastq.gz" --reference ref.fasta


Worfklows and subworkflows (--workflow | --subworkflow):

  - snippy                 Call Illumina PE variants with Snippy
      - snippy-core        Core genome SNP alignment with Snippy-Core
  - candidate              Candidate-guided SNP calls with Megalodon
  - denovo                 SNP calling with Clair or Medaka
  - random_forest          Random Forest SNP polishers
      - train              Train Random Forest classifiers on Snippy reference calls
      - evaluate           Evaluate Random Forests classifiers on sets of genomes

Help:

  nextflow run np-core/np-variants --workflow snippy --help true

Deployment and resource configuration:

  Specific process execution profiles defined within the configuration files are selected with
  the native argument `-profile`

  Resource configs can be configured hierarchically by first selecting a preset configuration
  file with `--config` or a specific resource preset configuration with `--resource_config`

  For more information see: https://github.com/np-core/config 

  Example:

      nextflow run np-core/np-variants --workflow snippy --config default -profile docker 

=========================================
```

## Core variants (Illumina)

`Snippy`

In its simplest incarnation, the variant calling workflow uses `Snippy` and `Gubbins` to generate a non-recombinant core-genome alignment for the provided set of isolates (includes output of non-core per-isolate variant calls) from high-quality short-read sequence data:

```
nextflow run np-core/np-variants --workflow core --fastq "isolates/*_R{1,2}.fq.gz" --reference ref.fasta
```

## Candidate Variants (ONT)

`Megalodon`

Candidate variants - for example core-genome SNPs called with `--workflow core` using `Snippy` - can be used as input to anchor the neural network basecalling output from `Guppy` and call candidates with `Megalodon`. `NanoPath` then merges candidates and calls into a filtered alignment which can be used for 'hybrid' phylogenetic trees, combining both Illumina and ONT sequence data. While this is useful to reconsruct divergences within a known evolutionary background, for example when producing barcoded nanopore panels from an outbreak of a known lineage for which sufficient background data is available, please note that within-outbreak branches would not consider novel variation (particularly if divergence is deep, that is, if the outbreak has been persisting and accumulated novel variation for some time) and thus, within-outbreak phylodynamic estimates based on candidate variants may be far from accurate.

`Guppy` is used for basecalling and should be run using GPU resourcing through configuration files and profiles - if you require more information regarding the selection of resources for this workflow, please visit: [`np-core/configs`](https://github.com/np-core/configs)

```
nextflow run np-core/np-variants --config jcu -profile tesla --workflow candidate --fast5 fast5/ --candidates core.vcf
```

## *De novo* Variants (ONT)

`Clair | Medaka`


## Random Forest SNP Polishers (ONT)

`Clair | Medaka`


SNP polishers trained with `Nextflow` sub-workflow: `train` 

The purpose of this workflow is to automate training random forest classifiers to remove false SNP calls from a VCF called with `Clair` or `Medaka`. We train the classifiers on referenced (Illumina PE) nanopore reads (barcoded, demultiplexed) over a range of coverage subsets (to account for coverage differences in the application of the classifiers).

In the preprint we show that performance is dependent on the reference genome chosen to train the classifiers, which should correspond to the reference genome you want to call the variants in your real data against. If high performance is desired, e.g. for functional SNP significance (removing false positives, antibiotic resistance related...) then the same reference genome as for the intended application, ideally from the same sequence type, should be selected and trained with this pipeline. Performance for phylogenetic applications can be more forgiving, and we have shown that even classifiers trained on a different species can remove sufficient false calls to allow accurate topological reconstruction and phylodynamic estimates. Performance of the same reference polisher critically depends on the training isolates (ONT + Illumina PE) - ideally they would belong to the same genotype / outbreak, but for all practical applications, performance is similarly high with isolates from the same sequence type, and even a combination of other sequence types of the same species.

We need the following ingredients for the classifiers:

* one or more reference sequences in `fasta` format to train against
* one or more collections of matching  ONT and Illumina PE reads from the same isolate

Training sets should be in `--train_dir` in named folders (which will be the prefixes for the trained models) with the following format:

```
train_dir/
  model_1/                   <-- model training collection
    isolate.fastq            <-- ONT
    isolate_R1.fastq.gz      <-- Illumina
    isolate_R2.fastq.gz      <-- Illumina
```

Execute the training pipeline on the training directory, using two reference genomes in the current directory (`jkd.fasta`, `dar.fasta`):

```
nextflow run np-core/np-variants -profile docker \
   --workflow random_forest \
   --subworkflow train \
   --caller clair \
   --train_dir training_sets \ 
   --outdir models \
   --train_references jkd.fasta,dar.fasta \
   --test_size 0.3
```

---

SNP polishers evaluated with `Nextflow` sub-workflow: `eval`

The purpose of this subworkflow is to evaluate the trained classifiers on sets of matching Illumina PE and ONT nanopore reads from the same isolates. This is mostly for the manuscript, but can be used to assess the performance of trained classifiers for their specific use cases if more matching short and long-read data on the same isolates is available.

We need the following ingredients for the evaluations:

* one or more reference sequences in `fasta` format to call reference and nanopore SNPs against
* one or more model files trained with the `train` subworkflow
* one or more collections of matching ONT and Illumina PE reads from the same isolate
