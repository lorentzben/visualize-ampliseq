process RENDERREPORT{

    tag "$reportName"
    label 'process_low'
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/microbiome_analyst:1.1' : 
    'lorentzb/microbiome_analyst:1.1' }"

    //container 'lorentzb/microbiome_analyst:1.1'

    input:

    val(reportName)
    path('results')
    path('metadata.tsv')
    path(report)
    path('item_of_interest.csv')
    path(table)

    output:

    path("*.html"), emit: html_report
    path("figures/*"), emit: figures_generated
     
    
    script:
    def table = table.name != 'NO_FILE' ? "$table" : ''
    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');
    Rscript -e "rmarkdown::render('01_report_MbA.Rmd', output_file='$PWD/01_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    #Rscript -e "rmarkdown::render('01_report_MbA.Rmd', output_file='$PWD/01_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}