#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hw2_pipeline } from './modules/local/hw2'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup'

params.reads = null
params.sra_id = null
params.samplesheet = null
params.reference = null
params.outdir = "results"
params.save_mpileup = false

process index_reference {
    tag "${meta.id}"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path(fasta), path("${fasta}.fai"), emit: reference

    script:
    """
    samtools faidx $fasta
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

process filter_variants {
    tag "${group_meta.group}"
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
        tuple val(group_meta), path(vcfs)

    output:
        tuple val(group_meta), path("${group_meta.group}_filtered"), emit: filtered

    script:
    def vcf_args = vcfs.join(' ')
    """
    mkdir ${group_meta.group}_filtered

    for vcf in ${vcf_args}; do
        name=\$(basename "\$vcf" .vcf.gz)

        bcftools view \\
            -i 'QUAL>=20' \\
            -Oz \\
            -o ${group_meta.group}_filtered/\${name}.filtered.vcf.gz \\
            "\$vcf"

        tabix -p vcf -f ${group_meta.group}_filtered/\${name}.filtered.vcf.gz
    done
    """

    stub:
    """
    mkdir ${group_meta.group}_filtered
    touch ${group_meta.group}_filtered/${group_meta.group}.filtered.vcf.gz
    touch ${group_meta.group}_filtered/${group_meta.group}.filtered.vcf.gz.tbi
    """
}

workflow {

    def input_count = [params.reads, params.sra_id, params.samplesheet].count { it }

    if (input_count != 1) {
        error "Please provide exactly one input: '--samplesheet', '--reads' or '--sra_id'"
    }

    if (!params.reference) {
        error "Please provide '--reference' to a genome fasta file"
    }

    hw2_pipeline()

    reference_for_index = Channel.value(
        tuple([id: 'reference'], file(params.reference, checkIfExists: true))
    )

    index_reference(reference_for_index)

    reference_input = index_reference.out.reference.first()

    bam_input = hw2_pipeline.out.bam.map { meta, bam, bai ->
        tuple(meta, bam, [], [])
    }

    BCFTOOLS_MPILEUP(bam_input, reference_input, params.save_mpileup)

    variants_by_group = BCFTOOLS_MPILEUP.out.vcf
        .map { meta, vcf ->
            tuple([group: meta.group], vcf)
        }
        .groupTuple()

    filter_variants(variants_by_group)
}