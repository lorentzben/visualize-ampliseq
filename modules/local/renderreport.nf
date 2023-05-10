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
    path('report.Rmd')
    file('item_of_interest.csv')
    path("table.tsv")
    path("table.qza")

    output:

    path("*.html"), emit: html_report
    path("figures/*"), emit: figures_generated
     
    
    script:
    
    '''
    #!/usr/bin/env bash

    ls -lRh
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');
    Rscript -e "rmarkdown::render('report.Rmd', output_file='$PWD/!{reportName}_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    #Rscript -e "rmarkdown::render('"!{report}"', output_file='$PWD/"!{report}"_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}