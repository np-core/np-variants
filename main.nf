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

        Model configuration files can be found inside the container at: /models

        Resources can be configured hierarchically by first selecting a configuration file from
        ${baseDir}/configs with `--config` and a specific resource configuration with `--resource config`

        Specific process execution profiles defined within the configuration files are selected with
        the native argument `-profile`

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

    Input / output configuration:

        --path                  the path to directory of fast5 files to pass to Guppy (single folder) or a glob for Fast5 [${params.path}]
        --archived              input files are expected to be tar gzipped files ending with .tar.gz or .tgz [${params.archived}]
        --outdir                the path to directory for result outputs [${params.outdir}]

    GPU basecalling configuration:

        --guppy_model               base guppy model configuration file for basecalling [${params.guppy_model}]
        --guppy_params              base guppy additional parameter configurations by user ["${params.guppy_params}"]
        --guppy_data                base guppy model data directory, inside container ["${params.guppy_data}"]
        --gpu_devices               gpus to use, provide list of devices passed to the -x flag in Guppy ["${params.gpu_devices}"]
        --gpu_forks                 parallel basecalling instances to launch on GPUs [${params.gpu_forks}]
        --runners_per_device        parallel basecalling runners on gpus [${params.runners_per_device}]
        --chunks_per_runner         the number of signal chunks processed on each runner [${params.chunks_per_runner}]
        --chunk_size                the size of the signal chunks processed on the gpu runers[${params.chunk_size}]
        --num_callers               the number of basecallers spread across the devices [${params.num_callers}]
        --cpu_threads_per_caller    the number of cpu threads per caller [${params.num_callers}]

    Qcat demultiplexing configuration:

        --demultiplex          activate demultiplexing with Qcat [${params.demultiplex}]
        --qcat_params          additional qcat parameters passed by the user ["${params.qcat_params}"]

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







// Declare external files
params.reference = "$PWD/ref.fasta"

params.meta_data = "$PWD/meta_data.tsv"

params.outdir = "$PWD/test_out"

include { Fastp } from './modules/fastp'
include { SnippyFastq } from './modules/snippy'
include { SnippyFasta } from './modules/snippy'
include { SnippyCore  } from './modules/snippy'
include { Gubbins  } from './modules/gubbins'

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

include { RAxML } from './modules/raxml'
include { TreeTime } from './modules/treetime'
include { VariantSites } from './modules/phybeast'
include { DateRandomisation } from './modules/phybeast'

// Basic phylogeny and phylodynamics based on ML

workflow phylodynamics_ml {
    take:
        alignment
    main:
        VariantSites(alignment)
        RAxML(VariantSites.out)
        TreeTime(RAxML.out, meta_data, alignment)
        DateRandomisation(RAxML.out, TreeTime.out[0], meta_data, alignment)
    emit:
        RAxML.out
}

// Advanced models on GPU using BEAST2 and BEAGLE

// include { BirthDeathSkyline } from './modules/beastling'
// include { CoalescentSkyline } from './modules/beastling'
// include { MultiTypeBirthDeath } from './modules/beastling'

// // Should be used for exploratory runs unless sure that priors are configured appropriately

// params.iterations = 1000000  
// params.coupled_mcmc = false
// params.beast_options = "--beagle_gpu"

// params.cosky_config = "default"
// params.cosky_dimensions = [2, 4, 8, 16]

// params.bdss_config = "default"
// params.bdss_dimensions = [5, 6, 7, 8]

// params.mtdb_config = "default"
// params.mtbd_types = ['mssa', 'clade']

// include { Beast as BeastCosky } from './modules/beast'
// include { Beast as BeastBDSS } from './modules/beast'
// include { Beast as BeastMTBD } from './modules/beast'

// workflow phylodynamics_beast {
//     take:
//         core_snp_alignment
//     main:
//         CoalescentSkyline() | BeastCosky
//         BirthDeathSkyline() | BeastBDSS
//         MultiTypeBirthDeath() | BeastMTBD
// }

// Execute

params.subworkflow = "megalodon"
params.panels = "$HOME/LINEAGES/ST93/Megalodon"
params.candidates = "$HOME/LINEAGES/ST93/core.vcf"
params.devices = "1"
params.guppy_server_path = "/opt-guppy/bin/guppy_basecall_server"
params.guppy_params = "--trim_barcodes --chunk_size 512 --chunks_per_runner 2048 --gpu_runners_per_device 4"
params.guppy_config = "dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg"

include { MegalodonVariants } from './modules/megalodon'

def get_barcode_fast5(dir){

    port = 5555
    return channel.fromPath("$params.panels/**/*", type: 'dir').map { port += 1; tuple(port, it.getParent().getName(), it.getName(), it) }
}

        
workflow megalodon_variants {
    take:
        barcodes
    main:
        MegalodonVariants(reference, barcodes)
    emit:
        MegalodonVariants.out
}

workflow {
    if (params.subworkflow == "assembly"){

        fasta = get_single_file("FNQ*.fasta") | snippy_fasta
        fastq = get_paired_files("*_{R1,R2}.fastq.gz") | snippy_fastq
        fasta.mix(fastq) | snippy_core | phylodynamics_ml
    
    } else if (params.subworkflow == "megalodon"){

        get_barcode_fast5(params.panels) | megalodon_variants

    }
}