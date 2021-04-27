#!/usr/bin/env nextflow


Channel.fromPath( params.in ).into { fastq_files }

println """\
         This workflow finds all the files matching your glob pattern and prints them. 
         Use it to check which files match to troubleshoot or prevent unexpected surprises. 
         Usage: nextflow run <path>/filepair_finder.nf --in="<glob>"
         Common glob patterns:
         "<path>/*.fastq.gz" for fastq.gz files
         Finding files in subfolders of a directory:
         "<path>/**/*.fastq.gz"
         File Finder for    NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}
         """
         .stripIndent()

fastq_files.subscribe { println it }

