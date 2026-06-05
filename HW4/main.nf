nextflow.enable.dsl = 2

include { HW2_PIPELINE } from './modules/local/hw2'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup'

params.accession = ''
params.samplesheet = ''
params.reference = ''
params.outdir = 'results'
params.threads = 4
params.save_mpileup = false

workflow {
    // Check inputs
    def input_count = [params.accession, params.samplesheet].count { it }
    if( input_count != 1 ) {
        error "Use exactly one input: --accession or --samplesheet"
    }
    if( !params.reference ) {
        error "Please provide --reference"
    }
    
    // Run pipeline
    HW2_PIPELINE()
    
    // Prepare BAM for variant calling
    bam_for_variants = HW2_PIPELINE.out.bam.map { meta, bam, bai, flagstat ->
        tuple(meta, bam, [], [])
    }
    
    // Prepare reference with index
    reference_with_index = HW2_PIPELINE.out.reference.map { meta, reference ->
        def fai = file("${reference}.fai")
        if (!fai.exists()) {
            "samtools faidx ${reference}".execute()
        }
        tuple(meta, reference, fai)
    }
    
    // Call variants
    BCFTOOLS_MPILEUP(bam_for_variants, reference_with_index, params.save_mpileup)
    
    // Group by group and filter
    variants_by_group = BCFTOOLS_MPILEUP.out.vcf
        .map { meta, vcf -> tuple([group: meta.group], vcf) }
        .groupTuple()
    
    FILTER_VARIANTS(variants_by_group)
}

process FILTER_VARIANTS {
    tag "${group_meta.group}"
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
    tuple val(group_meta), path(vcfs)

    output:
    tuple val(group_meta), path("${group_meta.group}_filtered")

    script:
    """
    mkdir ${group_meta.group}_filtered
    for vcf in ${vcfs}; do
        name=\$(basename "\$vcf" .vcf.gz)
        bcftools view -i 'QUAL>=20' -Oz -o ${group_meta.group}_filtered/\${name}.filtered.vcf.gz "\$vcf"
        tabix -p vcf -f ${group_meta.group}_filtered/\${name}.filtered.vcf.gz 2>/dev/null || true
    done
    """

    stub:
    """
    mkdir ${group_meta.group}_filtered
    touch ${group_meta.group}_filtered/dummy.vcf.gz
    """
}