process REPORT05ALPHABOXPLOT{
    
    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_05:2.0' : 
        'lorentzb/r_05:2.0' }"

    input: 

    val(reportName)
    path (vectors)
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('metadata.tsv')
    path('05_report.Rmd')

    output:

    path("05_report_*.html")
    path("alpha_diversity_boxplot/*")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('05_report.Rmd', output_file='$PWD/05_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    
    '''
}
