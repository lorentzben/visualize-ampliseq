process REPORT07RAREFACTION{

    tag "$reportName"
    label 'process_low'

    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_07:2.0' : 
        'lorentzb/r_07:2.0' }"

    input: 

    val(reportName)
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('07_report.Rmd')
    path(rarefact)
    path('metadata.tsv')

    output:

    path("07_report_*.html")
    path("07_report_*.pdf")
    path("rarefaction_plots/*")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}
