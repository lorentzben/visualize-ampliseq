process REPORT12PERMANOVA{
    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_12:2.0' : 
        'lorentzb/r_12:2.0' }"

    input:

    val(reportName)
    path('feature-table.qza')
    //TODO do we need this?
    path('results')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('taxonomy.qza')
    path('metadata.tsv')
    path(distances)
    path('12_report.Rmd')

    output:

    path("12_report_*.html"), emit: html_report
    path("anosim/*"), emit: anosim
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('12_report.Rmd', output_file='$PWD/12_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
    
    
    '''
}