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


// Workflow version

version = '0.1.1'

def helpMessage() {

    log.info"""
    =========================================
     N P - V A R I A N T S  v${version}
    =========================================

    Usage:

    A typical command for constructing the reference alignment and core-genome single nucleotide polymorphism calls:

        nextflow run np-core/np-variants --illumina fastq/ --tail "_R{1,2}.fastq.gz"

    Deployment and resource configuration:

        Resources can be configured hierarchically by:
            
            1. Selecting a tag or path to image file for Docker or Singularity (`--container`)
            2. Selecting a configuration file from ${baseDir}/configs (`--config`) 
            3. Selecting a resource configuration  ${baseDir}/configs (`--resource_config`)
            4. Selecting a native profile within the configuration file (`-profile`)

        --container             path to container file or docker tag to provision pipeline

                                  <np-core/signal>      Example for tag of Docker image
                                  <$HOME/phybeast.sif>  Example for path to Singularity image

        --config                select a configuration from the configs subdirectory of the pipeline

                                  <nextflow>  base configuration with docker or singularity profiles
                                  <jcu>       base configuration for the zodiac cluster at JCU
                                  <nectar>    base configuration for the nectar cluster at QCIF

        --resource_config       select a resource configuration nested within the selected configuration

                                  <process>   base resource configuration of processes for compute servers

        -profile                select a system executor profile from the config file - default:

                                  <docker> / <gpu_docker>  - expect container to be tag format
                                  <singularity> / <gpu_singularity> - expect container to be path to image

    Subworkflow selection:

        --workflow              select the variant subworkflow to select: core, candidate, denovo [${params.workflow}]
        --outdir                output directory for results from workflow [${params.outdir}]

    
    Subworkflow - Core Variants:

        --fastq | --fasta       glob to FASTQ and/or FASTA files for variant calling in Snippy ["${params.fastq}" | "${params.fasta}"]
        --reference             reference genome (FASTA) for alignment and variant calling [${params.reference}]

        --variant_sites         remove monomorphic sites from alignment after Gubbins [${params.variant_sites}]
        --snippy_params         string of additional parameters passed to Snippy ["${params.snippy_params}"]
        --gubbins_params        string of additional parameters passed to Gubbins ["${params.gubbins_params}"]

    Subworkflow - Megalodon:

        Model configuration files for Guppy can be found inside the container at: /models

        --fast5                 glob of directories containing Fast5 files for Megalodon [${prams.fast5}]
        --panels                path to nested panel directories, which contain barcode subdirectories (e.g. fast5/panel1/barcode01) [$params.panels]
        --candidates            VCF candidate variant file to select variants to call with Megalodon [${params.candidates}]

        --guppy_server_path     server path to guppy, inside container at: "/opt-guppy/bin/guppy_basecall_server" [${params.guppy_server_path}] 
        --guppy_params          string of additional parameters to pass to guppy, should contain model directory inside container: "-d /models" [${params.guppy_params}]
        --guppy_config          name of guppy basecalling configuration file [${params.guppy_config}]
        --guppy_devices         string of space delimited list of GPU devices for Guppy (e.g. "0 1") [${params.guppy_devices}]
    =========================================

    """.stripIndent()
}


params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def check_file(file) {
    
    path = Paths.get(file)

    if (path.exists()){
        log.info"""
        Detected input file: $file
        """
    } else {
        log.info"""
        Failed to detect input file: $file
        """
        exit 0
    }
}


if (params.recursive){
    _fastq_ont = ["${params.files}/**/*.fq", "${params.files}/**/*.fastq"]
    _fastq_illumina = ["${params.files}/**/*${params.tail}.fq.gz", "${params.files}/**/*${params.tail}.fastq.gz"]
} else {
    _fastq_ont = ["${params.files}/*.fq", "${params.files}/*.fastq"]
    _fastq_illumina = ["${params.files}/*${params.tail}.fq.gz", "${params.files}/*${params.tail}.fastq.gz"]
}

// Helper functions

def get_single_fastx(glob){
    return channel.fromPath(glob) | map { file -> tuple(file.baseName, file) }
}
def get_paired_fastq(glob){
    return channel.fromFilePairs(glob, flat: true)
}
def get_fast5(dir){
    return channel.fromPath("$params.fast5", type: 'dir').map { port += 1; tuple(port, it.getParent().getName(), it.getName(), it) }
}
def get_barcode_fast5(dir){
    port = 5555
    return channel.fromPath("$params.panels/**/*", type: 'dir').map { port += 1; tuple(port, it.getParent().getName(), it.getName(), it) }
}

// Workflow selection 
params.workflow = "core"


// Core (Illumina)

params.fastq = "*.fastq"
params.fasta = "*.fasta"

params.reference = "$PWD/ref.fasta"
check_file(params.reference)
reference = file(params.reference)  // stage the reference


// Candidates (Megalodon)
params.fast5 = ""
params.panels = "$HOME/LINEAGES/ST93/Megalodon"
params.candidates = "$HOME/LINEAGES/ST93/core.vcf"
params.devices = "1"
params.guppy_server_path = "/opt-guppy/bin/guppy_basecall_server"
params.guppy_params = "--trim_barcodes --chunk_size 512 --chunks_per_runner 2048 --gpu_runners_per_device 4"
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg"

params.outdir = "$PWD/test_out"

include { Fastp } from './modules/fastp'
include { SnippyFastq } from './modules/snippy'
include { SnippyFasta } from './modules/snippy'
include { SnippyCore  } from './modules/snippy'
include { Gubbins  } from './modules/gubbins'
include { Variants  } from './modules/PhyBeast'

workflow snippy_fastq {
    take:
        reads // id, forward, reverse
    main:
        Fastp(reads)
        SnippyFastq(Fastp.out, reference)
    emit:
        SnippyFastq.out // id, results
}       

workflow snippy_fasta {
    take:
        contigs // id, fasta
    main:
        SnippyFasta(contigs, reference)
    emit:
        SnippyFasta.out // id, results
}  

workflow snippy_core {
    take:
        snippy // results
    main:
        SnippyCore(snippy.collect(), reference)
        Gubbins(SnippyCore.out) // wgs snp alignment
    emit:
        Gubbins.out // non-recombinant snp core alignment
}



include { MegalodonVariants } from './modules/megalodon'



        
workflow megalodon_variants {
    take:
        barcodes
    main:
        MegalodonVariants(reference, barcodes)
    emit:
        MegalodonVariants.out
}

workflow {
    if (params.subworkflow == "core"){

        fasta = get_single_fastx(params.fasta) | snippy_fasta
        fastq = get_paired_fastq(params.fastq) | snippy_fastq
        fasta.mix(fastq) | snippy_core
    
    } else if (params.subworkflow == "candidate"){

        get_barcode_fast5(params.panels) | megalodon_variants

    }
}