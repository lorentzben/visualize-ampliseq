process REPORT04ALPHATABLE{

    tag "$reportName"
    label 'process_low'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_04:2.0' : 
        'lorentzb/r_04:2.0' }"


    input: 

    val(reportName)
    path (vectors)
    file('item_of_interest.csv')
    path('04_report.Rmd')

    output:

    path("04_report_*.html")
    path("04_report_*.pdf")
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    cp -L 04_report.Rmd $PWD/04_report_test.Rmd

    Rscript -e "rmarkdown::render('04_report_test.Rmd', output_file='04_report_$dt.html', output_format='html_document', output_dir='$PWD', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('04_report_test.Rmd', output_file='04_report_$dt.pdf', output_format='pdf_document', output_dir='$PWD', clean=TRUE, knit_root_dir='$PWD')"
    '''

}