process SRSCURVE{

    //tag "$table"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/srs:1.0' : 
        'lorentzb/srs:1.0' }"

    input:

    path('table.qza')
    path('table.tsv')
    path('srs_curve.rmd')
    path('my_count_table_min_max.py')
    

    output:

    path ("*.pdf"), emit: pdf_report
    path ("*.html"), emit: html_report
    path ("*.png"), emit: images 
    path("srs_max_curve_val.txt"), emit: max_val
    path("srs_min_curve_val.txt"), emit: min_val

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    '''
    #!/usr/bin/env bash

    python3 my_count_table_min_max.py
    
    Rscript -e "rmarkdown::render('srs_curve.rmd', output_file='$PWD/srs_curve.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    Rscript -e "rmarkdown::render('srs_curve.rmd', output_file='$PWD/srs_curve.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
    
    '''

}