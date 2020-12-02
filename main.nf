#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================================================================
                                        N P - S I G N A L   P I P E L I N E
==============================================================================================================================

 Nanopore signal processing pipeline (Fast5)

 Documentation: https://github.com/np-core/np-signal

 Original development by Queensland Genomics, Australian Institute of Tropical Health 
 and Medicince, The Peter Doherty Intitute for Infection and Immunity

Developers:

    Eike Steinig  @esteinig  < @EikeSteinig >

Pipeline part of the NanoPath core framework:

    https://github.com/np-core

NanoPath distributed pipeline framework Netflow:

    https://github.com/np-core/netflow

For interactive dashboard operation and report generation:

    https://github.com/np-core/nanopath

----------------------------------------------------------------------------------------
*/

import java.nio.file.Paths

nextflow.enable.dsl=2

// Helper functions

def check_path(p, descr) {
    
    if (p == ""){
        log.info"""
        No $descr provided (required)
        """
        exit 0
    } else {
        path = Paths.get(p)
        if (path.exists()){
            log.info"""
            Detected input $descr: $p
            """
        } else {
            log.info"""
            Failed to detect input $descr: $p
            """
            exit 0
        }
    }    
}


// Workflow selection 
params.workflow = "core"
params.outdir = "$PWD/results"

params.reference = "" // FASTA reference genome
if (params.reference){
    check_path(params.reference, "reference file") // required
    reference = file(params.reference)  // stage the reference
} else {
    reference = ""
}


// Core (Illumina)

params.fastq = "*_R{1,2}.fastq"
params.fasta = "" // optional
params.snippy_params = ""
params.gubbins_params = ""

// Candidates (Megalodon)

params.fast5 = ""
params.panels = ""

if (params.panels){
    check_path(params.panels, "panel directory")
}

params.candidates = "" // VCF
if (params.candidates){
    check_path(params.candidates, "candidate file")
    candidates = file(params.candidates) // stage file
} else {
    candidates = "" // define alt variable 
}

params.megalodon_params = ""
params.guppy_devices = "cuda:0"
params.reads_per_guppy_batch = 50

params.guppy_server_path = "/opt/ont/guppy/bin/guppy_basecall_server"  // reachable inside container
params.guppy_params = "-d /guppy_models" // should always include "-d /guppy_models" or "-d /rerio_models/"
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg" // Rerio: res_dna_r941_min_modbases-all-context_v001.cfg

// De novo

params.caller = "clair"
params.medaka_model = "r941_min_high_g360"
params.clair_model = "/clair_models/model"
params.clair_haploid = "--haploid_sensitive"
params.coverage = ""
params.genome_size = "2.8m"

// Forest Evaluation

params.dir_snippy = ""
params.dir_ont = ""
params.eval_mask_weak = 0.8
params.eval_models = ""
params.eval_caller = "clair"

if ( params.eval_models && params.eval_models instanceof String ){
    eval_models = params.eval_models.split(",").collect { file(it) }
} else {
    eval_models = params.eval_models
}

if ( params.coverage instanceof String ){
    coverage = params.coverage.split(",").collect { it }
} else {
    coverage = params.coverage
}

// Workflow version

version = '0.1.4'

def helpMessage() {

    log.info"""
    =========================================
     N P - V A R I A N T S  v${version}
    =========================================

    Usage:

    A typical command to construct a core-genome variant alignment from PE Illumina data:

        nextflow run np-core/np-variants --workflow core --path "*_R{1,2}.fastq.gz" --reference ref.fasta

    Deployment and resource configuration:

            Resources can be configured hierarchically by first selecting a configuration file from
            presets with `--config` and a resource presets with `--resource_config`

            Specific process execution profiles defined within the configuration files are selected with
            the native argument `-profile`

            For more information see: https://github.com/np-core/config 

    Available subworkflows:

        - snippy | snippy-core     call variants from Illumina PE with Snippy / core genome alignment with Snippy-Core
        - candidate_megalodon
        - candidate_clair
        - denovo
        - forest_training
        - forest_evaluation
        - forest_polish

    Required configuration:

        --workflow                  select the variant subworkflow to select: core, candidate, denovo [${params.workflow}]        
        --reference                 reference genome (FASTA) for alignment and variant calling [${params.reference}]
        --outdir                    output directory for results from workflow [${params.outdir}]

    Snippy / Snippy-Core:

        --fastq | --fasta           glob to FASTQ and/or FASTA files for variant calling in Snippy ["${params.fastq}" | "${params.fasta}"]
        --snippy_params             string of additional parameters passed to Snippy ["${params.snippy_params}"]
        --gubbins_params            string of additional parameters passed to Gubbins ["${params.gubbins_params}"]

    Candidate Megalodon:

        Model configuration files for Guppy can be found inside the container at: /guppy_models

        --fast5                     directory or glob of directories containing Fast5 files for a single sample [${params.path}]
        --panels                    path to nested panel directories, which contain barcode subdirectories (e.g. fast5/panel1/barcode01) [$params.panels]
        --candidates                VCF candidate variant file to select variants to call with Megalodon [${params.candidates}]

        --reads_per_guppy_batch     number of reads to batch for concurrent processing with Guppy [${params.reads_per_guppy_batch}]
        --guppy_server_path         server path to guppy, inside container at: "/usr/bin/local" [${params.guppy_server_path}] 
        --guppy_params              string of additional parameters to pass to Guppy, should contain model directory inside container: "-d /guppy_models" [${params.guppy_params}]
        --guppy_config              name of guppy basecalling configuration file [${params.guppy_config}]
        --guppy_devices             string of space delimited list of GPU devices for Guppy (e.g. "0 1") [${params.guppy_devices}]

    DeNovo Variants

        --fastq                     glob to FASTQ files for variant calling with Medaka or Clair ["${params.fastq}"]
        --caller                    variant caller to use, one of: medaka, clair [${params.caller}]
        --medaka_model              Medaka model to use for variant calling [${params.medaka_model}]
        --coverage                  Comma delimited string of target coverage to subsample before variant calling [${params.subsample}]
        --genome_size               Genome size for subsampling in Rasusa [${params.genome_size}]

    SNP Polisher: Random Forest Training

        --
    
    SNP Polisher: Random Forest Evaluation

        --

    SNP Polisher: Random Forest Polishing

        --

    =========================================

    """.stripIndent()
}


params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Helper functions

def get_single_file(glob){
    channel.fromPath(glob, type: 'file') | map { path -> tuple(path.baseName, path) }
}
def get_paired_fastq(glob){
    return channel.fromFilePairs(glob, flat: true)
}
def get_fast5_dir(dir){
    return channel.fromPath(dir, type: 'dir').map { tuple(it.getName(), it) }
}
def get_fast5_panels(dir){
    return channel.fromPath("${dir}/**/*", type: 'dir').map { tuple(it.getParent().getName(), it.getName(), it) }
}
def get_evaluation_batches(snippy_dir, ont_dir){

    ont = Channel.fromFilePairs("${ont_dir}/*.{vcf,txt}", flat: true, type: 'file')
    
    snippy = Channel.fromPath("${snippy_dir}/*.ref.vcf").map { tuple(it.simpleName, it) }
    
    data = ont.join(snippy) | map { tuple(it[0], it[3], it[2], it[1]) }

    return data

}

include { Fastp } from './modules/fastp'
include { SnippyFastq } from './modules/snippy'
include { SnippyFasta } from './modules/snippy'
include { SnippyCore  } from './modules/snippy'
include { Gubbins  } from './modules/gubbins'
include { MegalodonVariants } from './modules/megalodon'
include { MegalodonVariantsPanels } from './modules/megalodon'
include { MedakaVariants } from './modules/medaka'
include { MinimapONT } from './modules/minimap2'
include { ClairVariants } from './modules/clair'
include { RasusaMulti } from './modules/rasusa'

include { TrainRandomForest } from './modules/variants'
include { RasusaMultiTraining } from './modules/rasusa'
include { MinimapMultiTraining } from './modules/minimap2'
include { ClairVariantsTraining } from './modules/clair'
include { MedakaVariantsTraining } from './modules/medaka'

include { EvaluateRandomForest } from './modules/variants'
include { ProcessRandomForestEvaluations } from './modules/variants'

workflow snippy_fastq {
    take:
        reads // id, forward, reverse
    main:
        Fastp(reads)
        SnippyFastq(Fastp.out, reference)
    emit:
        SnippyFastq.out // results
}       

workflow snippy_fasta {
    take:
        contigs // id, fasta
    main:
        SnippyFasta(contigs, reference)
    emit:
        SnippyFasta.out // results
}  

workflow snippy_core {
    take:
        snippy // results
    main:
        SnippyCore(snippy.collect(), reference)
        Gubbins(SnippyCore.out[0]) // wgs snp alignment
    emit:
        Gubbins.out // non-recombinant snp core alignment
}

        
workflow megalodon_dir {
    take:
        fast5  // id, f5d
    main:
        MegalodonVariants(reference, candidates, fast5)
    emit:
        MegalodonVariants.out
}

workflow megalodon_panels {
    take:
        panels // panel, barcode, f5d
    main:
        MegalodonVariantsPanels(reference, candidates, panels)
    emit:
        MegalodonVariantsPanels.out
}

workflow denovo_snps {
    take:
        fastq // id, fq
    main:
        if (params.coverage){
            fastq = RasusaMulti(fastq, coverage)
        }
        mapped = MinimapONT(fastq, reference)
        if (params.caller == "medaka"){
            variants = MedakaVariants(mapped, reference)
        } else if (params.caller == "clair"){
            variants = ClairVariants(mapped, reference)
        }
    emit:
        variants[0]  // vcf variant files
        variants[1] // bam alignments
}


workflow evaluate_forest {
    take:
        eval_batch  // per file: id, snippy_vcf, ont_vcf, ont_stats
    main:
        EvaluateRandomForest(eval_batch, eval_models) | collect | ProcessRandomForestEvaluations
    emit:
        null
}

include { TrainingReferenceSnippy } from './modules/snippy'

params.training_set_illumina_glob = "*_R{1,2}.fastq.gz"
params.training_set_ont_glob = "*.fastq"

params.train_dir = ""
params.test_size = 0.3

params.train_references = ""
train_references = params.train_references.split(",").collect { it }

train_coverages = [2, 5, 10, 30, 50, 100]

def showTrainingConfiguration() {
    log.info"""

    Model Training
    ==============

    Variant caller    : ${params.caller}
    Model directory   : ${params.train_dir}
    References        : ${train_references}
    Coverage subsets  : ${train_coverages}
    

    Model evaluation
    ================


    """.stripIndent()
} 


def get_train_data(train_dir){

    // Get the training data from train_dir/{model_name}

    illumina = channel.fromFilePairs("${train_dir}/**/${params.training_set_illumina_glob}", type: 'file', flat: true)
        .map { tuple(it[1].getParent().getName(), it[0],  it[1], it[2]) } // model, id, fw, rev

    ont = channel.fromPath("${train_dir}/**/${params.training_set_ont_glob}", type: 'file')
        .map { tuple(it.getParent().getName(), it.simpleName,  it) } // model, id, fq

    return illumina.join(ont, by: [0, 1])

}

workflow train_forest {
    take:
        train_data  // model_name, isolate_id, reference_name, reference_file, ont_fq, illumina_vcf
    main:
        fastq_model_cov = RasusaMultiTraining(train_data, train_coverages)
        mapped_model_cov = MinimapMultiTraining(fastq_model_cov)
        if (params.caller == "medaka"){
            variants_model_cov = MedakaVariantsTraining(mapped_model_cov)
        } else if (params.caller == "clair"){
            variants_model_cov = ClairVariantsTraining(mapped_model_cov)
        }
        variants_model_cov | groupTuple(by: [0, 1] ) | TrainRandomForest   // by model_name, reference        
    emit:
        null

}

workflow publication {

    take:
        reads
    main:
        // Step 1: Training isolate sets: call Illumina reference VCFs for each reference genome (in param.dir_train: set1/*.ref.fastq, set1/*.ont.fastq)

        get_train_data(params.train_dir) | view


        // Step2: Training isolate sets: train RF polishers for each reference genome


        // Step 3: Evaluation isolate sets: call the evaluation isolate reference Illumina VCFs with Snippy for each reference genome



        // Step 4:Evaluation isolate sets: evaluate the different RF models on ONT data in comparison to Illumina reference VCFs

    emit:
        null

}

workflow {
    
    if (params.workflow == "candidate_megalodon"){
        if (params.panels){
            println "Calling Megalodon candidate SNPs ($params.caller) on multiplex Fast5 panel input directory: $params.panels"
            get_fast5_panels(params.panels) | megalodon_panels
        } else {
            println "Calling Megalodon candidate SNPs ($params.caller) on Fast5 input directory: $params.panels"
            get_fast5_dir(params.fast5) | megalodon_dir
        }

    } else if (params.workflow == "snippy" | params.workflow == "snippy-core") {
        println "Calling SNPs using Snippy on input: $params.fastq (fastq) / $params.fasta (fasta)"
        if (params.fastq){
            fastq = get_paired_fastq(params.fastq) | snippy_fastq
        } else {
            fastq = channel.empty()
        }    
        if (params.fasta){
            fasta = get_single_file(params.fasta) | snippy_fasta
        } else {
            fasta = channel.empty()
        }

        snippy = fasta.mix(fastq)
        
        if (params.workflow == "snippy-core"){
            println "Calling core SNPs using SnippyCore and Gubbins"
            snippy | snippy_core
        }
    
    } else if (params.workflow == "denovo"){
        println "Calling de novo SNPs ($params.caller) on FASTQ input: $params.fastq"
        get_single_file(params.fastq) | denovo_snps
    } else if (params.workflow == "forest_evaluation"){
        println "Forest evaluation parsed reference variant calls in directory: $params.dir_snippy"
        println "Forest evaluation parsed ONT variant calls in directory: $params.dir_ont"
        get_evaluation_batches(params.dir_snippy, params.dir_ont) | evaluate_forest
    } else if (params.workflow == "forest_training"){
        println "Forest training using collections in directory: $params.train_dir"
        get_train_data(params.train_dir) | train_forest
    } else if (params.workflow == "publication"){
        showTrainingConfiguration()
        train_data = get_train_data(params.train_dir)
        TrainingReferenceSnippy(train_data, train_references)

    }


}