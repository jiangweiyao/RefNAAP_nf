# Reference Based Nanopore Assembly Analysis Pipeline (RefNAAP) for Nextflow

This pipeline processes the raw, basecalled Nanopore reads generated from the Nanopore Sequencing. It does the following steps:
1. It QCs the fastq files using fastqc and multiqc to generate a quality report.
2. It trims the left and right ends of the reads by 25 basepairs, and filters out reads shorter than 50bp.
4. It performs reference based assembly of the reads, fixes gaps, and polishes using Medaka.

# How to run
Nextflow version > 20 is requires
To run locally using singularity:
```
nextflow run 'https://github.com/jiangweiyao/RefNAAP_nf' -profile singularity 
```

A dummy awsbatch profile is include in the repo. Link it to your AWS Batch environment with your credentials to run it.

## Author
- Jiangwei Yao
- Crystal Gigante

## License 
Apache License 2.0

