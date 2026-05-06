#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hw2_pipeline } from './modules/local/hw2'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup'

params.reads = null
params.sra_id = null
params.reference = null
params.outdir = "results"

if (params.reads && params.sra_id) {
    error "Please provide either '--reads' or '--sra_id', not both"
}
if (!params.reads && !params.sra_id) {
    error "Please provide either '--reads' or '--sra_id'"
}
if (!params.reference) {
    error "Please provide '--reference' to a genome fasta file"
}

workflow {
    // 1. Run the main pipeline (from HW2)
    hw2_pipeline()
    
    // 2. Prepare input data for bcftools mpileup
    //    Format: tuple(meta, bam, intervals_mpileup, intervals_call)
    bam_input = hw2_pipeline.out.bam.map { sample, bam, bai ->
        tuple([id: sample], bam, [], [])
    }
    
    // 3. Reference with index
    //    Create .fai if it doesn't exist
    ref_fai = file(params.reference + ".fai")
    if (!ref_fai.exists()) {
        "samtools faidx ${params.reference}".execute()
    }
    
    reference_input = hw2_pipeline.out.reference.map { ref ->
        tuple([id: ref.baseName], ref, ref_fai)
    }
    
    // 4. Call variants using nf-core module
    BCFTOOLS_MPILEUP(bam_input, reference_input, false)
}