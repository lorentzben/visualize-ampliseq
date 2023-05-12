process REPORT08RANKEDABUNDANCE{

    tag "$reportName"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_08:2.0' : 
        'lorentzb/r_08:2.0' }"

    input: 
     
    val(reportName)
    path('feature-table.qza')
    path('rooted-tree.qza')
    path('taxonomy.qza')
    path('metadata.tsv')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('08_report.Rmd')
    
    
    output:

    path("08_report_*.html")
    path("08_report_*.pdf")
    path("ranked_abundance_curves/*")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}
