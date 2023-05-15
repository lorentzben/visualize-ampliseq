process REPORT09UNIFRACHEATMAP{

    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_09:2.0' : 
        'lorentzb/r_09:2.0' }"

    input:

    val(reportName)
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path("metadata.tsv")
    path(distances)
    path("09_report.Rmd")

    output:

    path("09_report_*.html"), emit: html_report
    path("09_report_*.pdf"), emit: pdf_report
    path("unifrac_heatmaps/*"), emit: unifrac_png
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}