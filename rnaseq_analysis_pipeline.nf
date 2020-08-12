#!/usr/bin/env nextflow

log.info """
-----------------------------------
- Rnaseq pipeline analysis        -
- version 0.1                     -
- Created by: Adrian Odrzywolski  -
-----------------------------------



"""
def help_message(String append_message = null){
    log.info """
    ${append_message}

    How to use pipeline:
    nextflow run rnaseq_pipeline_analysis.nf --reads '*_R{1,2}.fastq.gz' --genome '~/reference/GRCh37.fasta' --with-docker
    
    Required arguments:
    --reads - path to single- or paired-end reads in form of regular exporession. File can be gziped (perferebly)
    --genome - path to reference genome
    --annotation - path to annotation

    Optional arguments:
    --output - select directory path where all output needs to be stored [default: .]
    """.stripIndent()
}

// help message
if(params.help){
    help_message()
    exit 0
}

if(!params.reads || !params.genome){
    help_message("Reads or reference genome path missing!")
    exit 1
}

// making basic channels - these are crucial!

    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { exit 1, help_message("No pairs of reads were found: ${params.reads}")}
        .into { reads_raw_fastqc; reads_trim_raw_reads }

/* 
    TODO: single end!
    Channel
        .fromPath(params.reads as String, checkIfExists: true)
        .ifEmpty { exit 1, help_message("No singles of reads were found: ${params.reads}")}
        .into { reads_raw_fastqc; reads_trim_raw_reads } */


Channel
    .fromPath(params.genome as String, checkIfExists: true)
    .ifEmpty { exit 1, help_message("No reference was found: ${params.genome}")}
    .set { genome }

Channel
    .fromPath(params.annotation as String, checkIfExists: true)
    .ifEmpty { exit 1, help_message("No annotation file was found: ${params.annotation}")}
    .into {annotation_bq ; annotation_hisat2}


// Assuming that there is no index build

process reference_index {
    input:
    file reference_path from genome
   
    output:
    file("${reference_path.baseName}*ht2") into hisat2_index

    script:
    
    """
    hisat2-build -p ${task.cpus} -f ${reference_path} ${reference_path.baseName} 
    """
    
}

// Doing FASTQC on raw reads

process raw_fastqc {
    cpus = 2
    publishDir params.output, mode: 'copy', overwrite: false, pattern: '*.{zip,html}'

    input:
    tuple sampleId, file(reads_path) from reads_raw_fastqc

    output:
    file "*.{zip,html}" into result_raw_fastqc

    script:
    """
    fastqc $reads_path -t ${task.cpus} --noextract --nogroup --outdir .
    """
}

// Trimming bbduk

process trim_raw_reads {
    input:
    tuple sampleId, file(reads_path) from reads_trim_raw_reads
    output:
    tuple sampleId, file("*_trimmed*") into result_trim_raw_reads, result_trim_raw_reads_align

    script:
    """
    reads_1=${reads_path[0]}
    reads_1_simple_name=${reads_path[0].simpleName}
    reads_1_extension="\${reads_1#*.}"
    reads_2=${reads_path[1]}
    reads_2_simple_name=${reads_path[1].simpleName}
    reads_2_extension="\${reads_2#*.}"

    bash /bio-tools/bbmap/bbduk.sh -Xmx5g threads=${task.cpus} \\
    in=\$reads_1 \\
    in2=\$reads_2 \\
    out=\$reads_1_simple_name"_trimmed."\$reads_1_extension \\
    out2=\$reads_2_simple_name"_trimmed."\$reads_2_extension \\
    trimq=6 \\
    maq=10 \\
    ref=/bio-tools/bbmap/resources/adapters.fa \\
    k=23 \\
    mink=11 \\
    hdist=1 \\
    ktrim=w \\
    tbo=t \\
    tpe=t \\
    overwrite=t \\
    stats=bbduk_info.txt 2> stats_trimming.txt
    """
}

// Doing FASTQC on trimmed reads

process trimmed_fastqc {
    cpus = 2
    
    publishDir params.output, mode: 'copy', overwrite: false, pattern: '*.{zip,html}'
    input:
    tuple sampleId, file(reads_path) from result_trim_raw_reads

    output:
    file "*.{zip,html}" into result_trimmed_fastqc

    script:
    """
    fastqc $reads_path -t ${task.cpus} --noextract --nogroup --outdir .
    """
}

// Doing aligment with hisat2

process hisat2{
    publishDir params.output, mode: 'copy', overwrite: false, pattern: '*.txt'
    input:
    file reference from hisat2_index.collect()
    tuple sampleId, file(reads_path) from result_trim_raw_reads_align

    output:
    file "*bam" into fresh_bams
    file "*_summary.txt" into hisat2_summary

    script:
    """
    reads_1=${reads_path[0]}
    reads_1_simple_name=${reads_path[0].simpleName}
    reads_1_extension="\${reads_1#*.}"
    reads_2=${reads_path[1]}
    reads_2_simple_name=${reads_path[1].simpleName}
    reads_2_extension="\${reads_2#*.}"

    hisat2 \\
    -p ${task.cpus} \\
    --summary-file ${sampleId}_Hisat_summary.txt \\
    -x ${reference[0].simpleName} \\
    -1 \${reads_1} \\
    -2 \${reads_2} \\
    | samtools view -hbS > ${sampleId}.bam """
}

// Fixing mate with samtools

process fixmate_bam {
    input:
    file bam from fresh_bams

    output:
    file "*_fixmate.bam" into fixed_bam

    script:
    """
    samtools fixmate -@ ${task.cpus} -O bam ${bam} ${bam.baseName}_fixmate.bam
    """
}

// Sorting bam with samtools

process sort_bam {
    input:
    file bam from fixed_bam

    output:
    file "*_fixmate_sorted.bam" into fixed_sorted_bam, fixed_sorted_bam_fc,fixed_sorted_bam_bamqual,fixed_sorted_bam_raport

    script:
    """
    samtools sort ${bam} -@ ${task.cpus} -o ${bam.baseName}_sorted.bam
    """
}

// Making index of bam

process index_bam {
    cpus = 1
    input:
    file bam from fixed_sorted_bam

    output:
    file "*_fixmate_sorted.bam.bai" into indexof_fixed_sorted_bam

    script:
    """
    samtools index ${bam}
    """
}

// Counting features with subread

process feature_count {
    publishDir params.output, mode: 'copy', overwrite: false
    
    input:
    file annotation_path from annotation_hisat2
    file bam from fixed_sorted_bam_fc
    file bam_index from indexof_fixed_sorted_bam

    output:
    file '*.summary' into feature_count_result

    script:
    """
    featureCounts -p -B -T ${task.cpus} -t exon \\
    -a ${annotation_path} \\
    -o ${bam.simpleName} \\
    ${bam}
    """
}

// Making different quality matricies with qualimap

process bam_qual {
    publishDir params.output, mode: 'copy', overwrite: true
    input:
    file annotation_path from annotation_bq
    file bam from fixed_sorted_bam_bamqual
    output:
    file "${bam.simpleName}_stats/*" into bam_qual_result
    script:
    """
    qualimap --java-mem-size=8G -outdir ${bam.simpleName}_stats --bam ${bam} -pe -s -gtf ${annotation_path}
    """
}

// Final raport with multiqc

process make_report {
    publishDir params.output, mode: 'copy', overwrite: true

    input:
    file bam from fixed_sorted_bam_raport
    file bam_qual from bam_qual_result
    file hisat2 from hisat2_summary
    file raw_fastqc from result_raw_fastqc
    file trimmed_fastqc from result_trimmed_fastqc
    file feature_count from feature_count_result

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"


    script:
    """
    multiqc -f -c /bio-tools/config.yaml .
    """
}


workflow.onComplete {
    if(workflow.success){
        log.info """
        Done!
        """
        multiqc_report
    }else{
        log.info """
        Run completed with errors!
        """
    }
}

