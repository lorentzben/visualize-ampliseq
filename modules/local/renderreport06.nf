
process REPORT06ORDINATION{

    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_06:2.0' : 
        'lorentzb/r_06:2.0' }"

    input: 

    val(reportName)
    path('feature-table.qza')
    path('rooted-tree.qza')
    path('taxonomy.qza')
    path('metadata.tsv')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('06_report.Rmd')
    path(pcoas)
    path(vectors)

    output:

    path("06_report_*.html")
    path("06_report_*.pdf")
    path("beta_diversity_ordination/*")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REPORT06BNMDSORDINATION{

    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_06:2.0' : 
        'lorentzb/r_06:2.0' }"

    input: 

    val(reportName)
    path('feature-table.qza')
    path('rooted-tree.qza')
    path('taxonomy.qza')
    path('metadata.tsv')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('06b_report.Rmd')
    path(pcoas)
    path(vectors)

    output:

    path("06b_report_*.html")
    path("06b_report_*.pdf")
    path("beta_diversity_ordination/*")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('06b_report.Rmd', output_file='$PWD/06b_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('06b_report.Rmd', output_file='$PWD/06b_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}
