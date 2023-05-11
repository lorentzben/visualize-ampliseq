process REPORT03HEATMAP{

    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_03:2.0' : 
        'lorentzb/r_03:2.0' }"

    input: 

    val(reportName)
    path('feature-table.qza')
    path('rooted-tree.qza')
    path('taxonomy.qza')
    path('metadata.tsv')
    path('03_report.Rmd')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')
    path('overall_summary.tsv')

    output:

    path("03_report_*.html"), emit: report_html
    path("03_report_*.pdf"), emit: report_pdf
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    cp -L 03_report.Rmd $PWD/03_report_test.Rmd

    Rscript -e "rmarkdown::render('03_report_test.Rmd', output_file='$PWD/03_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('03_report_test.Rmd', output_file='$PWD/03_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}