process REPORT10BETABOXPLOT{
    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_10:2.0' : 
        'lorentzb/r_10:2.0' }"

    input:

    val(reportName)
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path("metadata.tsv")
    path("10_report.Rmd")
    path(distances)

    output:

    path("10_report_*.html"), emit: html_report
    path("10_report_*.pdf"), emit: pdf_report
    path("beta_div_boxplots/*"), emit: beta_boxplots
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''

}