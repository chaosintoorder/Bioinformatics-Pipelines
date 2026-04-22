params.reads = null
params.sra_id = null
params.reference = null
params.outdir = "results"

// =======================
// PROCESS: FETCH READS
// =======================
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

// =======================
// PROCESS: FASTQC
// =======================
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

// =======================
// PROCESS: TRIMMING
// =======================
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

// =======================
// WORKFLOW: TRIMMED QC
// =======================
workflow trimmed_qc_wf {
    take:
        trimmed_reads
    main:
        fastqc('trimmed', trimmed_reads)
    emit:
        fastqc.out
}

// =======================
// PROCESS: ASSEMBLY
// =======================
process assembly {
    tag "$sample_id"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    spades.py -t ${task.cpus} -1 $r1 -2 $r2 -o spades_out
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}

// =======================
// PROCESS: MAPPING
// =======================
process mapping {
    tag "$sample_id"
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
        tuple val(sample_id), path(r1), path(r2), path(ref)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa index $ref
    bwa mem -t ${task.cpus} $ref $r1 $r2 | \
        samtools view -@ ${task.cpus} -bS - | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    """
}

// =======================
// PROCESS: COVERAGE PLOT
// =======================
process plot_coverage {
    tag "$sample_id"
    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam)
        path r_script  // принимаем R скрипт как входной файл
    
    output:
        tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")
    
    script:
    """
    samtools depth $bam > ${sample_id}_depth.txt
    Rscript $r_script ${sample_id}
    """
}

// =======================
// MAIN WORKFLOW
// =======================
workflow {
    if (params.sra_id) {
        raw_reads_ch = fetch_reads(params.sra_id)
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Please provide either '--sra_id' or '--reads \"data/*_{1,2}.fastq\"'"
    }
    
    fastqc('raw', raw_reads_ch)
    
    trimmed_reads_ch = trimm(raw_reads_ch)
    
    trimmed_qc_wf(trimmed_reads_ch)
    
    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference, checkIfExists: true)

        mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    } else {
        assembly_ch = assembly(trimmed_reads_ch)
        mapping_in_ch = trimmed_reads_ch.join(assembly_ch)
    }
    
    bam_ch = mapping(mapping_in_ch)
    
    plot_coverage(bam_ch, file("plot_coverage.R"))
}