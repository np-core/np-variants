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

def check_path(p) {
    
    path = Paths.get(p)

    if (path.exists()){
        log.info"""
        Detected input path: $p
        """
    } else {
        log.info"""
        Failed to detect input path: $p
        """
        exit 0
    }
}


// Workflow selection 
params.workflow = "core"
params.outdir = "$PWD/results"


params.reference = "" // FASTA reference genome
check_path(params.reference) // required
reference = file(params.reference)  // stage the reference


// Core (Illumina)
params.fastq = "*_R{1,2}.fastq"
params.fasta = "" // optional
params.snippy_params = ""
params.gubbins_params = ""

// Candidates (Megalodon)

params.path = ""
params.panels = ""

if (params.panels){
    check_path(params.panels)
}

params.candidates = "ref.fasta" // VCF
if (params.candidates){
    check_path(params.candidates)
    candidates = file(params.candidates) // stage file
} else {
    candidates = "" // define alt variable 
}

params.megalodon_params = ""
params.gpu_devices = "1"
params.reads_per_guppy_batch = 50

params.guppy_server_path = "/opt/ont/guppy/bin/guppy_basecall_server"  // reachable inside container
params.guppy_params = "-d /guppy_models" // should always include "-d /guppy_models" or "-d /rerio_models/" or "/.../barcoding" models
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg" // Rerio: res_dna_r941_min_modbases-all-context_v001.cfg



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
        --outdir                    output directory for results from workflow [${params.outdir}]

    Subworkflow - Core Variants:

        --fastq | --fasta           glob to FASTQ and/or FASTA files for variant calling in Snippy ["${params.fastq}" | "${params.fasta}"]
        --reference                 reference genome (FASTA) for alignment and variant calling [${params.reference}]

        --variant_sites             remove monomorphic sites from alignment after Gubbins [${params.variant_sites}]
        --snippy_params             string of additional parameters passed to Snippy ["${params.snippy_params}"]
        --gubbins_params            string of additional parameters passed to Gubbins ["${params.gubbins_params}"]

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

    =========================================

    """.stripIndent()
}


params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Helper functions

def get_single_fasta(glob){
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

include { Fastp } from './modules/fastp'
include { SnippyFastq } from './modules/snippy'
include { SnippyFasta } from './modules/snippy'
include { SnippyCore  } from './modules/snippy'
include { Gubbins  } from './modules/gubbins'
include { MegalodonVariants } from './modules/megalodon'
include { MegalodonVariantsPanels } from './modules/megalodon'

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

workflow {
    
    if (params.workflow == "candidate"){
        // ONT candidate workflow with Megalodon
        if (params.panels){
            get_fast5_panels(params.panels) | megalodon_panels
        } else {
            get_fast5_dir(params.path) | megalodon_dir
        }

    } else if (params.workflow == "core") {
        // Illumina core-genome SNP workflow with Snippy and Gubbins
        if (params.fastq){
            fastq = get_paired_fastq(params.fastq) | snippy_fastq
        } else {
            fastq = channel.empty()
        }    
        if (params.fasta){
            fasta = get_single_fasta(params.fasta) | snippy_fasta
        } else {
            fasta = channel.empty()
        }
    
        fasta.mix(fastq) | snippy_core
    
    } else if (params.workflow == "denovo"){



    }


}