process BCFTOOLS_MPILEUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(intervals_mpileup), path(intervals_call)
    tuple val(meta2), path(fasta), path(fai)
    val save_mpileup

    output:
    tuple val(meta), path("*vcf.gz"), emit: vcf
    tuple val(meta), path("*vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*stats.txt"), emit: stats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --output-type u \\
        ${bam} \\
        | bcftools call --multiallelic-caller --variants-only \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z

    tabix -p vcf ${prefix}.vcf.gz
    bcftools stats ${prefix}.vcf.gz > ${prefix}.stats.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.stats.txt
    """
}
