#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "${projectDir}/ampliseq_results"
params.metadata = "${projectDir}/metadata.tsv"
params.ioi = "treatment"
params.outdir = "results"

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input }
         metadata : ${params.metadata}
         item of interest : ${params.ioi}
         outdir   : ${params.outdir}
         """
         .stripIndent()


process REPORT01BARPLOT{

    publishDir "results/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "results/html", pattern: "*.html", mode: "copy"

    container 'docker://lorentzb/r_01:2.0'


    input:

    path core-metric
    path rooted-tree
    path taxonomy
    path metadata
    path order-item-of-interest
    path 01-report
    val item-of-interest
    

    script:

    '''
    #!/usr/bin/env bash

    echo !{item-of-interest} > item_of_interest.csv
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('01_report.Rmd', output_file='$PWD/01_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('01_report.Rmd', output_file='$PWD/01_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''

}

workflow{
    input_ch = Channel.fromPath(params.input)
    metadata_ch = Channel.fromPath(params.metadata)
    core-metric = Channel.fromPath(params.input+'/qiime2/')
    REPORT01BARPLOT(input_ch)
}