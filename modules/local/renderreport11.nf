process REPORT11UPGMA{

    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_11:2.0' : 
        'lorentzb/r_11:2.0' }"

    input:

    val(reportName)
    path('feature-table.qza')
    path('rooted-tree.qza')
    path('taxonomy.qza')
    path('metadata.tsv')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('11_report.Rmd')

    output:

    path("11_report_*.html"), emit: html_report
    path("11_report_*.pdf"), emit: pdf_report
    path("upgma_plots/*"), emit: upgma_plots
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}