# HW4: Multi-sample NGS Pipeline with CSV Input and Variant Filtering

==Deadline: 05.06==
Performed by Victor Florin

Start with the pipeline created in HW3, and add the following things:  
1. Add support for multi-samples input via .csv file (use this and previous  
lecture materials) (several viruses or several batteries as groups of  
samples)  
1.1. Read all samples to one channel  
1.2. Split them by some field  
1.3. Process (to call variants)  
1.4. Join back to one channel  
2. Run final analysis on new process – e.g. filter variants  
3. For the new process use stub running for designing and implementing.


I need to create a folder called HW4 and copy the contents of HW3 into it.
I deleted the 'work' folder in HW3, so I will redownload the reads for SRR1175163.

Since a multi-sample input is expected, I have chosen SRR1175163 (HW2, HW3), SRR38739167 and SRR38739168. As these are all **dengue virus** paired-end samples, I'll group them by experimental design: SRR1175163 = randomPCR, SRR38739167 and SRR38739168 = targetCapture.

I copied the contents of the HW3 folder into a new folder called HW4, created a reads folder inside HW4, and then downloaded the SRA datasets:
```
$ prefetch SRR1175163 -v
$ prefetch SRR38739167 -v
$ prefetch SRR38739168 -v
```

Next, I obtained the .fastq files:
```
$ fasterq-dump --split-files SRR1175163/SRR1175163.sra
$ fasterq-dump --split-files SRR38739167/SRR38739167.sra
$ fasterq-dump --split-files SRR38739168/SRR38739168.sra
```

Removed the .sra files:
```
$ rm -rf SRR1175163
$ rm -rf SRR38739167
$ rm -rf SRR38739168
```

As the reference genome I used the [Dengue virus 1 complete genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1?report=fasta), downloaded it and placed it into the HW4 folder under the name reference.fasta.

Created the samplesheet.csv:
```
sample,group,fastq_1,fastq_2
SRR1175163,randomPCR,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR1175163_1.fastq,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR1175163_2.fastq
SRR38739167,targetCapture,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR38739167_1.fastq,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR38739167_2.fastq
SRR38739168,targetCapture,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR38739168_1.fastq,/home/student/Bioinformatics-Pipelines/HW4/reads/SRR38739168_2.fastq
```

Modified the modules/local/hw2.nf file:
```
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
```

Adjusted the contents of modules/nf-core/bcftools/mpileup/main.nf:
```
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
```

Edited the root main.nf:
```
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
```

Updated nextflow.config:
```
// =======================
// DEFAULT PARAMETERS
// =======================
params {
    reads = null
    sra_id = null
    samplesheet = null
    reference = null
    outdir = "results"
    save_mpileup = false
    conda_env = null
    container = 'chaosintoorder/hw3-pipeline:latest'
}

// =======================
// DEFAULT PROCESS SETTINGS
// =======================
process {
    cpus = 4
    memory = '8 GB'
    time = '4h'
    errorStrategy = 'retry'
    maxRetries = 3

    withName: 'BCFTOOLS_MPILEUP' {
        publishDir = [
            path: { "${params.outdir}/variants" },
            mode: 'copy'
        ]
    }

    withName: 'index_reference' {
        publishDir = [
            path: { "${params.outdir}/reference" },
            mode: 'copy'
        ]
    }

    withName: 'filter_variants' {
        publishDir = [
            path: { "${params.outdir}/filtered_variants" },
            mode: 'copy'
        ]
    }

    withName: 'plot_coverage' {
        publishDir = [
            path: { "${params.outdir}/coverage" },
            mode: 'copy'
        ]
    }

    withName: 'fastqc' {
        publishDir = [
            path: { "${params.outdir}/qc" },
            mode: 'copy'
        ]
    }

    withName: 'fetch_reads' {
        publishDir = [
            path: { "${params.outdir}/raw_reads" },
            mode: 'copy'
        ]
    }

    withName: 'trimm' {
        publishDir = [
            path: { "${params.outdir}/trimmed_reads" },
            mode: 'copy',
            pattern: "*_p.fastq"
        ]
    }
}

// =======================
// PROFILES
// =======================
profiles {
    local {
        process {
            cpus = 4
            memory = '2 GB'
            time = '8h'
            conda = params.conda_env ?: '/home/student/miniconda3/envs/hw3'
        }
        conda.enabled = true
        docker.enabled = false
    }

    cluster {
        process {
            executor = 'slurm'
            cpus = 16
            memory = '64 GB'
            time = '48h'
            queue = 'normal'
            conda = "${projectDir}/environment.yml"
        }
        conda.enabled = true
        docker.enabled = false
    }

    container {
        process {
            cpus = 4
            memory = '16 GB'
            time = '8h'
            container = params.container
        }
        conda.enabled = false
        docker.enabled = true
        docker.runOptions = '--rm --platform linux/amd64'
    }
}
```

Added htslib to the root environment.yml:
```
name: hw3-pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc=0.12.1
  - trimmomatic=0.39
  - sra-tools=3.0.0
  - spades=3.15.5
  - bwa=0.7.17
  - samtools=1.16.1
  - bcftools=1.17
  - htslib=1.17
  - r-base=4.2
  - r-ggplot2
```

Checked the versions of the required tools:
```
nextflow -version
fastqc --version
trimmomatic -version
bwa 2>&1 | head
samtools --version | head
bcftools --version | head
Rscript --version
```

Trying the following command:
```
nextflow run . \
-profile local \
--conda_env /home/student/miniconda3/envs/hw3 \
--samplesheet samplesheet.csv \
--reference /absolute/path/to/reference.fasta \
--outdir results_stub \
-stub-run
```