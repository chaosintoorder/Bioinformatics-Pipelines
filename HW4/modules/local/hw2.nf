process fetch_reads {
    tag "${meta.id}"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    
    input:
        val meta

    output:
        tuple val(meta), path("${meta.id}_1.fastq"), path("${meta.id}_2.fastq")

    script:
    """
    fasterq-dump --split-files ${meta.id}
    """
}

process fastqc {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        val reads_type
        tuple val(meta), path(r1), path(r2)

    output:
        path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}

process trimm {
    tag "${meta.id}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq"
    
    input:
        tuple val(meta), path(r1), path(r2)

    output:
        tuple val(meta), path("${meta.id}_R1_p.fastq"), path("${meta.id}_R2_p.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} \\
        $r1 $r2 \\
        ${meta.id}_R1_p.fastq ${meta.id}_R1_u.fastq \\
        ${meta.id}_R2_p.fastq ${meta.id}_R2_u.fastq \\
        LEADING:3 TRAILING:3 MINLEN:36
    """
}

workflow trimmed_qc_wf {
    take:
        trimmed_reads

    main:
        fastqc('trimmed', trimmed_reads)

    emit:
        fastqc.out
}

process alignment {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
        tuple val(meta), path(r1), path(r2), path(ref)

    output:
        tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai")

    script:
    """
    bwa index $ref

    bwa mem -t ${task.cpus} $ref $r1 $r2 | \\
        samtools view -@ ${task.cpus} -bS - | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    """
}

process plot_coverage {
    tag "${meta.id}"
    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
        tuple val(meta), path(bam), path(bai)
        path r_script
    
    output:
        tuple val(meta), path("${meta.id}_coverage.png"), path("${meta.id}_depth.txt")
    
    script:
    """
    samtools depth $bam > ${meta.id}_depth.txt
    Rscript $r_script ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_coverage.png
    touch ${meta.id}_depth.txt
    """
}

workflow hw2_pipeline {
    main:
    if (params.samplesheet) {
        raw_reads_ch = Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    [id: row.sample, group: row.group],
                    file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true)
                )
            }
    } else if (params.sra_id) {
        raw_reads_ch = fetch_reads(Channel.of([id: params.sra_id, group: params.sra_id]))
    } else if (params.reads) {
        raw_reads_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .map { sample_id, reads ->
                tuple(
                    [id: sample_id.toString(), group: sample_id.toString()],
                    reads[0],
                    reads[1]
                )
            }
    } else {
        error "Please provide one input: '--samplesheet', '--sra_id' or '--reads'"
    }
    
    fastqc('raw', raw_reads_ch)
    
    trimmed_reads_ch = trimm(raw_reads_ch)
    
    trimmed_qc_wf(trimmed_reads_ch)
    
    ref_ch = Channel.value(file(params.reference, checkIfExists: true))
    
    mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    bam_ch = alignment(mapping_in_ch)
    
    plot_coverage(bam_ch, file("${projectDir}/bin/plot_coverage.R"))
    
    emit:
    bam = bam_ch
    reference = ref_ch
}