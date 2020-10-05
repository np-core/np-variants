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

params.path = ""
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
params.guppy_params = "-d /guppy_models" // should always include "-d /guppy_models" or "-d /rerio_models/" or "/.../barcoding" models
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg" // Rerio: res_dna_r941_min_modbases-all-context_v001.cfg

// De novo

params.caller = "medaka"
params.medaka_model = "r941_min_high_g360"
params.clair_model = "/clair_models/model"
params.clair_haploid = "--haploid_sensitive"
params.coverage = ""
params.genome_size = "2.8m"

// Forest Evaluation

params.dir_snippy = ""
params.dir_ont = ""
params.eval_mask_weak = 0.8
params.eval_models = "$baseDir/random_forest/test_model.composite.sav"
params.eval_caller = "clair"

if ( params.eval_models instanceof String ){
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

    Subworkflow selection:

        --workflow                  select the variant subworkflow to select: core, candidate, denovo [${params.workflow}]        
        --reference                 reference genome (FASTA) for alignment and variant calling [${params.reference}]
        --outdir                    output directory for results from workflow [${params.outdir}]

    Subworkflow - Core Variants:

        --fastq | --fasta           glob to FASTQ and/or FASTA files for variant calling in Snippy ["${params.fastq}" | "${params.fasta}"]
        --snippy_params             string of additional parameters passed to Snippy ["${params.snippy_params}"]
        --gubbins_params            string of additional parameters passed to Gubbins ["${params.gubbins_params}"]

    Subworkflow - DeNovo Variants

        --fastq                     glob to FASTQ files for variant calling with Medaka or Clair ["${params.fastq}"]
        --caller                    variant caller to use, one of: medaka, clair [${params.caller}]
        --medaka_model              Medaka model to use for variant calling [${params.medaka_model}]
        --coverage                  Comma delimited string of target coverage to subsample before variant calling [${params.subsample}]
        --genome_size               Genome size for subsampling in Rasusa [${params.genome_size}]

    Subworkflow - Megalodon Haploid Variants:

        Model configuration files for Guppy can be found inside the container at: /guppy_models

        --path                      directory or glob of directories containing Fast5 files for a single sample [${params.path}]
        --panels                    path to nested panel directories, which contain barcode subdirectories (e.g. fast5/panel1/barcode01) [$params.panels]
        --candidates                VCF candidate variant file to select variants to call with Megalodon [${params.candidates}]

        --reads_per_guppy_batch     number of reads to batch for concurrent processing with Guppy [${params.reads_per_guppy_batch}]
        --guppy_server_path         server path to guppy, inside container at: "/usr/bin/local" [${params.guppy_server_path}] 
        --guppy_params              string of additional parameters to pass to Guppy, should contain model directory inside container: "-d /guppy_models" [${params.guppy_params}]
        --guppy_config              name of guppy basecalling configuration file [${params.guppy_config}]
        --guppy_devices             string of space delimited list of GPU devices for Guppy (e.g. "0 1") [${params.guppy_devices}]

    Subworkflow - Random Forest SNP Evaluation

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
def get_fast5_panel(dir){
    return channel.fromPath("${dir}/**/*", type: 'dir').map { tuple(it.getParent().getName(), it.getName(), it) }
}
def get_evaluation_batches(snippy_dir, ont_dir){

    snippy_vcf = Channel.fromPath("${snippy_dir}/*.vcf", type: 'file').map { tuple(it.baseName, it) }
    
    ont_vcf = Channel.fromPath("${ont_dir}/*.vcf", type: 'file').map { tuple(it.baseName, it) }
    ont_stats = Channel.fromPath("${ont_dir}/*.txt", type: 'file').map { tuple(it.baseName, it) }

    ont = ont_vcf.cross(ont_stats).unique().map { crossed ->
        if (crossed[0][0] == crossed[1][0]){ // id same
            return tuple( crossed[0][0], crossed[0][1], crossed[1][1] )  // id, ont_vcf, stats
        } 
    }

    matches = snippy_vcf.cross(ont).map { crossed ->
        if (crossed[0][0] == crossed[1][0]){ // id same
            return tuple( crossed[0][0], crossed[0][1], crossed[1][1], crossed[1][2] )   // id, snippy_vcf, ont_vcf, stats
        } 
    }

    return matches

}
def get_train_collections(snippy_dir, ont_dir){

    snippy_vcf = Channel.fromPath("${snippy_dir}/*.vcf", type: 'file').map { tuple(it.baseName, it) }
    
    ont = Channel.fromFilePairs("${ont_dir}/**/*.{vcf,txt}", type: 'file', flat: true).map { tuple(it[0], it[1].getParent().getName(),  it[1], it[2]) }
    
    matches = snippy_vcf.mix(ont) | groupTuple
    
    matches | view
    
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
include { EvaluateRandomForest } from './modules/variants'
include { ProcessRandomForestEvaluations } from './modules/variants'
include { TrainRandomForest } from './modules/variants'

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

workflow train_forest {
    take:
        train_data  // per file: train_collection_id, snippy_vcf, ont_vcf, ont_stats
    main:
        TrainRandomForest(train_data)
    emit:
        null

}

workflow evaluate_forest {
    take:
        eval_batch  // per file: id, snippy_vcf, ont_vcf, ont_stats
    main:
        EvaluateRandomForest(eval_batch, eval_models) | collect | ProcessRandomForestEvaluations
    emit:
        null

}

workflow {
    
    if (params.workflow == "candidate"){
        // ONT candidate workflow with Megalodon
        if (params.panels){
            get_fast5_panels(params.panels) | megalodon_panels
        } else {
            get_fast5_dir(params.path) | megalodon_dir
        }

    } else if (params.workflow == "core") {
        // Illumina core-genome SNP workflow with FASTQ / FASTA --> Snippy, SnippyCore and Gubbins
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
    
        fasta.mix(fastq) | snippy_core
    
    } else if (params.workflow == "denovo"){
        // ONT denovo workflow with Medaka or Clair
        get_single_file(params.fastq) | denovo_snps
    } else if (params.workflow == "forest_evaluation"){
        // Random Forest Classifier Evaluation
        get_evaluation_batches(params.dir_snippy, params.dir_ont) | evaluate_forest
    } else if (params.workflow == "forest_training"){
        println "Forest training"
        get_train_collections(params.dir_snippy, params.dir_ont)
    }


}