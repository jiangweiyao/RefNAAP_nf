#!/usr/bin/env nextflow

//params.in = "$baseDir/fastq/*.fastq.gz"
//params.out = "$HOME/eefnaap_nf_out"
//params.reference = "$baseDir/reference/Americas2.fasta"
//params.model = "r10_min_high_g303"
//params.size = 50
//params.left = 25
//params.right = 25

//AmericansN = file("$baseDir/reference/AmericasN.bed")
//AmericansG = file("$baseDir/reference/AmericasG.bed")

reference = Channel.fromPath(params.reference).into{ reference1; reference2 }
model = params.model
fastq_files = Channel.fromPath(params.in, type: 'file').map { file -> tuple(file.simpleName, file) }
fastq_files.into { fastq_files1; fastq_files2 }

scaffold_cutter_ch = Channel.fromPath("$baseDir/scaffold_cutter.R")
gapfixer_ch = Channel.fromPath("$baseDir/gapfixer.R")


println """\
         RabiesRefNAAP     NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}
         outdir (--out)                      : ${params.out}
         reference file (--reference)        : ${params.reference}
         medaka model (--model)              : ${params.model}
         Minimum read length (--size)        : ${params.size}
         Trim left n bases (--left)          : ${params.left}
         Trim right n bases (--right)        : ${params.right}
         """
         .stripIndent()


process seqtk_trim_filter {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    memory '3 GB'

    input:
    tuple name, file(fastq) from fastq_files2

    output:
    tuple name, file("*_filtered.fq") into fastq_filtered

    """
    seqtk trimfq -b ${params.left} -e ${params.right} $fastq > ${name}_trimmed.fq
    seqtk seq -L ${params.size} ${name}_trimmed.fq > ${name}_filtered.fq
    """

}

process fastqc {
    
    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    memory '3 GB' 


    input:
    set val(name), file(fastq) from fastq_files1 
 
    output:
    file "*_fastqc.{zip,html}" into qc_files
    file "*_fastqc.{zip,html}" into qc_files1

    """
    fastqc ${fastq}
    """
}


process multiqc {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    memory '3 GB'

    input:
    file reports  from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into final1

    """
    multiqc .
    """
}

process minimap2 {

    errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple val(name), file(fastq) from fastq_filtered
    file ref from reference1.first()

    output:
    tuple file(fastq), file("${fastq.simpleName}.sam") into sam_files

    """
    minimap2 -ax map-ont $ref ${fastq} > ${fastq.simpleName}.sam
    """
}

process samtools {

    errorStrategy 'ignore'
    publishDir params.out, pattern: "*_cov.txt", mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple file(fastq), file(sam_file)from sam_files

    output:
    tuple file(fastq), path("${sam_file.simpleName}.coverage") into coverage_files
    //path("${sam_file.simpleName}_cov.txt") into coverage_stat 

    """
    samtools view -S -b ${sam_file.simpleName}.sam > ${sam_file.simpleName}.bam
    samtools sort ${sam_file.simpleName}.bam -o ${sam_file.simpleName}.sorted.bam
    samtools index ${sam_file.simpleName}.sorted.bam
    samtools depth ${sam_file.simpleName}.sorted.bam > ${sam_file.simpleName}.coverage

    """
}

process scaffold {

    errorStrategy 'ignore'
    //publishDir params.out, pattern: "${coverage_file.simpleName}_cov.txt", mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple file(fastq), file(coverage_file) from coverage_files
    file ref from reference2.first()
    file(scaffold_cutter) from scaffold_cutter_ch.first()

    output:
    tuple file(fastq), path("${coverage_file.simpleName}_scaffold.fasta") into scaffold_files
    path "${coverage_file.simpleName}_cov.txt" into cov_stat_files

    """
    Rscript $scaffold_cutter 1 10 $ref ${coverage_file.simpleName}.coverage ${coverage_file.simpleName}_scaffold.fasta ${coverage_file.simpleName}_cov.txt   
    """
}

process medaka {

    errorStrategy 'ignore'
//    publishDir params.out, mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple file(fastq), file(scaffold_file) from scaffold_files

    output:
    tuple file(fastq), file(scaffold_file), path("${fastq.simpleName}_medaka/consensus.fasta") into medaka_files

    """
    medaka_consensus -i $fastq -d $scaffold_file -o ${fastq.simpleName}_medaka -m $model 
    """
}

process gapfixer {

    errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple file(fastq), file(scaffold_file), file(medaka_file) from medaka_files
    file(gapfixer) from gapfixer_ch.first()

    output:
    tuple file(fastq), path("${scaffold_file.simpleName}_gapfixed.fasta") into final_files

    """
    Rscript $gapfixer $scaffold_file $medaka_file ${scaffold_file.simpleName}_gapfixed.fasta
    """
}


process consensus_call {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    memory '7 GB'

    input:
    tuple file(fastq), file(scaffold_file) from final_files

    output:
    path("${fastq.simpleName}_final.fasta") into consensus_output
    //path("${fastq.simpleName}_medaka3") into consensus_output
    """
    medaka_consensus -i $fastq -d $scaffold_file -o ${fastq.simpleName}_medaka2 -m $model
    bcftools mpileup -f $scaffold_file ${fastq.simpleName}_medaka2/calls_to_draft.bam | bcftools call -mv -Oz -o ${fastq.simpleName}.vcf.gz --ploidy 1
    bcftools index ${fastq.simpleName}.vcf.gz
    cat $scaffold_file | bcftools consensus ${fastq.simpleName}.vcf.gz -H 1 > ${fastq.simpleName}_int.fasta
    medaka_consensus -i $fastq -d ${fastq.simpleName}_int.fasta -o ${fastq.simpleName}_medaka3 -m $model
    cp ${fastq.simpleName}_medaka3/consensus.fasta ${fastq.simpleName}_final.fasta
    """
}
