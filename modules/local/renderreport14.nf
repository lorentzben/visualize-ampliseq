process REPORT14CITATIONS{

    tag "$reportName"
    label 'process_low'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_14:2.0' : 
        'lorentzb/r_14:2.0' }"

    input:
    val(reportName)
    path("14_report.Rmd")

    output:

    path("14_report_*.html"), emit: html_report
    path("14_report_*.pdf"), emit: pdf_report
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''


}