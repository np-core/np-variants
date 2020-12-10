#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================================================================
                                        N P - V A R I A N T S   P I P E L I N E
==============================================================================================================================

 Nanopore variant calling and processing pipeline (Fast5)

 Documentation: https://github.com/np-core/np-variants

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
import groovy.io.*

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
params.subworkflow = ""
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
params.clair_model = "/clair_models/model"
params.clair_haploid = "--haploid_sensitive"
params.coverage = ""
params.genome_size = "2.8m"


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



// PUBLICATION REPLICATION

// Global

params.illumina_glob = "*_R{1,2}.fastq.gz"
params.ont_glob = "*.fastq"

// Training 
params.train_dir = ""
params.test_size = 0.3
params.train_references = ""
if (params.train_references) {
    train_references = params.train_references.split(",").collect { file(it) }
} else {
    train_references = null
}

train_coverages = [2, 5, 10, 30, 50, 100]

params.eval_dir = ""
params.model_dir = ""
params.mask_weak = 0.8
params.eval_references = ""

if (params.eval_references) {
    eval_references = params.eval_references.split(",").collect { file(it) }
} else {
    eval_references = null
}

// Evaluation

model_collections = ""
evaluation_collections = ""

include { FastpTraining } from './modules/fastp'
include { SnippyTraining } from './modules/snippy'
include { RasusaTraining } from './modules/rasusa'
include { MinimapTraining } from './modules/minimap2'
include { ClairTraining } from './modules/clair'
include { MedakaTraining } from './modules/medaka'
include { RandomForestTraining } from './modules/variants'

include { FastpEvaluation } from './modules/fastp'
include { SnippyEvaluation  } from './modules/snippy'
include { MinimapEvaluation  } from './modules/minimap2'
include { ClairEvaluation } from './modules/clair'
include { MedakaEvaluation } from './modules/medaka'
include { EvaluateRandomForest } from './modules/variants'
include { ProcessEvaluations } from './modules/variants'



def showTrainingConfiguration() {
    log.info"""
    
    Global config
    ==============

    Variant caller          : ${params.caller}
    Oxford Nanopore         : ${params.ont_glob}
    Illumina                : ${params.illumina_glob}

    Model training
    ==============

    Model directory         : ${params.train_dir}
    Training references     : ${train_references}
    Coverage subsets        : ${train_coverages}
    Training test size      : ${params.test_size}
    
    ---------------
    
    Model training sets     : ${model_collections}

    """.stripIndent()
} 

def showEvaluationConfiguration(){
    
    log.info"""
    
    Global config
    ==============

    Variant caller          : ${params.caller}
    Oxford Nanopore         : ${params.ont_glob}
    Illumina                : ${params.illumina_glob}

    Model evaluation
    ================

    Model directory        : ${params.model_dir}
    Evaluation directory   : ${params.eval_dir}
    Evaluation references  : ${params.eval_references}
    Mask weak sites        : ${params.mask_weak}

    """.stripIndent()

}

def get_train_data(train_dir){

    // Get the training data from train_dir/{model_name}

    illumina = channel.fromFilePairs("${train_dir}/**/${params.illumina_glob}", type: 'file', flat: true)
        .map { tuple(it[1].getParent().getName(), it[0], it[1], it[2]) } // model, id, fw, rev

    ont = channel.fromPath("${train_dir}/**/${params.ont_glob}", type: 'file')
        .map { tuple(it.getParent().getName(), it.simpleName,  it) } // model, id, fq

    return illumina.join(ont, by: [0, 1])

}

def get_eval_illumina(eval_dir){

    // Get the training data from train_dir/{model_name}

    return channel.fromFilePairs("${eval_dir}/**/${params.illumina_glob}", type: 'file', flat: true)
        .map { tuple(it[1].getParent().getName(), it[0], it[1], it[2]) } // eval set, id, fw, rev

}

def get_eval_ont(eval_dir){

    // Get the training data from train_dir/{model_name}

    return channel.fromPath("${eval_dir}/**/${params.ont_glob}", type: 'file')
        .map { tuple(it.getParent().getName(), it.simpleName,  it) } // eval set, id, fq

}


params.medaka_model = "r941_min_high_g360"

if ( file(params.medaka_model).exists() ){
    params.medaka_model  = file(params.medaka_model)
}


workflow train_forest {
    // model_name, isolate_id, reference_name, reference_file, ont_fq, illumina_vcf

    train_data = get_train_data(params.train_dir) | FastpTraining
    models = SnippyTraining(train_data, train_references) 

    fastq_model_cov = RasusaTraining(models, train_coverages)
    mapped_model_cov = MinimapTraining(fastq_model_cov)

    if (params.caller == "medaka"){
        variants_model_cov = MedakaTraining(mapped_model_cov)
    } else if (params.caller == "clair"){
        variants_model_cov = ClairTraining(mapped_model_cov)
    }

    variants_model_cov | groupTuple(by: [0, 1] ) | RandomForestTraining   // by model_name, refname, reference     
}


workflow eval_forest {
   
    illumina = get_eval_illumina(params.eval_dir) | FastpEvaluation
    illumina_snps = SnippyEvaluation(illumina, eval_references)  // call reference Illumina evaluation isolates with Snippy for each reference

    ont = get_eval_ont(params.eval_dir)
    mapped = MinimapEvaluation(ont, eval_references) // ONT SNP calls with Clair for each reference

    if (params.caller == "medaka"){
        ont_snps = MedakaEvaluation(mapped)
    } else if (params.caller == "clair"){
        ont_snps = ClairEvaluation(mapped)
    }

    snps = illumina_snps.join(ont_snps, by: [0, 1, 2])  // merge by (id, evaluation set, reference) matches (ont / illumina)

    eval_models = Channel.fromPath("${params.model_dir}/*.composite.sav", type: 'file')
    EvaluateRandomForest(snps, eval_models) | collect | ProcessEvaluations

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

    } else if (params.workflow == "snippy" | params.workflow == "snippy_core") {
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
        
        if (params.workflow == "snippy_core"){
            println "Calling core SNPs using SnippyCore and Gubbins"
            snippy | snippy_core
        }
    
    } else if (params.workflow == "denovo"){
        println "Calling de novo SNPs ($params.caller) on FASTQ input: $params.fastq"
        get_single_file(params.fastq) | denovo_snps
    } else if (params.workflow == "random_forest"){

        if (params.subworkflow == "train") {
            showTrainingConfiguration() 
            train_forest()
        } else if (params.subworkflow == "eval"){
            showEvaluationConfiguration()
            eval_forest()
        }

    }


}

// """

//     Required configuration:

//         --workflow                  select the variant subworkflow to select: core, candidate, denovo [${params.workflow}]        
//         --reference                 reference genome (FASTA) for alignment and variant calling [${params.reference}]
//         --outdir                    output directory for results from workflow [${params.outdir}]

//     Snippy / Snippy-Core:

//         --fastq | --fasta           glob to FASTQ and/or FASTA files for variant calling in Snippy ["${params.fastq}" | "${params.fasta}"]
//         --snippy_params             string of additional parameters passed to Snippy ["${params.snippy_params}"]
//         --gubbins_params            string of additional parameters passed to Gubbins ["${params.gubbins_params}"]

//     Candidate Megalodon:

//         Model configuration files for Guppy can be found inside the container at: /guppy_models

//         --fast5                     directory or glob of directories containing Fast5 files for a single sample [${params.path}]
//         --panels                    path to nested panel directories, which contain barcode subdirectories (e.g. fast5/panel1/barcode01) [$params.panels]
//         --candidates                VCF candidate variant file to select variants to call with Megalodon [${params.candidates}]

//         --reads_per_guppy_batch     number of reads to batch for concurrent processing with Guppy [${params.reads_per_guppy_batch}]
//         --guppy_server_path         server path to guppy, inside container at: "/usr/bin/local" [${params.guppy_server_path}] 
//         --guppy_params              string of additional parameters to pass to Guppy, should contain model directory inside container: "-d /guppy_models" [${params.guppy_params}]
//         --guppy_config              name of guppy basecalling configuration file [${params.guppy_config}]
//         --guppy_devices             string of space delimited list of GPU devices for Guppy (e.g. "0 1") [${params.guppy_devices}]

//     DeNovo Variants

//         --fastq                     glob to FASTQ files for variant calling with Medaka or Clair ["${params.fastq}"]
//         --caller                    variant caller to use, one of: medaka, clair [${params.caller}]
//         --medaka_model              Medaka model to use for variant calling [${params.medaka_model}]
//         --coverage                  Comma delimited string of target coverage to subsample before variant calling [${params.subsample}]
//         --genome_size               Genome size for subsampling in Rasusa [${params.genome_size}]

//     SNP Polisher: Random Forest Training

//         --
    
//     SNP Polisher: Random Forest Evaluation

//         --

//     SNP Polisher: Random Forest Polishing

//         --

// """