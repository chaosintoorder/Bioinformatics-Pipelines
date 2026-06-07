process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(bam), path(intervals_mpileup), path(intervals_call)
    tuple val(meta2), path(fasta), path(fai)
    val save_mpileup
    
    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${meta.id}.stats.txt"), emit: stats

    script:
    """
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --output-type u \\
        ${bam} \\
        | bcftools call --multiallelic-caller --variants-only \\
        | bcftools view --output-file ${meta.id}.vcf.gz --output-type z
    
    tabix -p vcf -f ${meta.id}.vcf.gz
    bcftools stats ${meta.id}.vcf.gz > ${meta.id}.stats.txt
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.vcf.gz
    touch ${meta.id}.vcf.gz.tbi
    touch ${meta.id}.stats.txt
    """
}