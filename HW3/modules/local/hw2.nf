process fetch_reads {
    tag "$sra_id"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    
    input:
        val sra_id

    output:
        tuple val(sra_id), path("${sra_id}_1.fastq"), path("${sra_id}_2.fastq")

    script:
    """
    fasterq-dump --split-files $sra_id
    """
}

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/qc_${reads_type}", mode: 'copy'
    
    input:
        val reads_type
        tuple val(sample_id), path(r1), path(r2)

    output:
        path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}

process trimm {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq"
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_R1_p.fastq"), path("${sample_id}_R2_p.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        $r1 $r2 \
        ${sample_id}_R1_p.fastq ${sample_id}_R1_u.fastq \
        ${sample_id}_R2_p.fastq ${sample_id}_R2_u.fastq \
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
    tag "$sample_id"
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2), path(ref)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa index $ref
    bwa mem -t ${task.cpus} $ref $r1 $r2 | \
        samtools view -@ ${task.cpus} -bS - | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process plot_coverage {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam), path(bai)
        path r_script
    
    output:
        tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")
    
    script:
    """
    samtools depth $bam > ${sample_id}_depth.txt
    Rscript $r_script ${sample_id}
    """
}

workflow hw2_pipeline {
    main:
    if (params.sra_id) {
        raw_reads_ch = fetch_reads(params.sra_id)
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Please provide either '--sra_id' or '--reads'"
    }
    
    // raw QC
    fastqc('raw', raw_reads_ch)
    
    // trimming
    trimmed_reads_ch = trimm(raw_reads_ch)
    
    // trimmed QC (через workflow обертку)
    trimmed_qc_wf(trimmed_reads_ch)
    
    // reference
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    // mapping
    mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    bam_ch = alignment(mapping_in_ch)
    
    // coverage
    plot_coverage(bam_ch, file("${projectDir}/bin/plot_coverage.R"))
    
    emit:
    bam = bam_ch
    reference = ref_ch
}