process REPORT13LEFSE{
    tag "$reportName"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_13:2.0' : 
        'lorentzb/r_13:2.0' }"

    input:

    val(reportName)
    path("lefse_images/*")
    path('13_report.Rmd')
    path('13_report_local.Rmd')
    file('item_of_interest.csv')
    path('order_item_of_interest.csv')

    output:

    path("13_report_*.html"), emit: html_report
    
         

    script:
    if(workflow.profile.contains('local'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('13_report_local.Rmd', output_file='$PWD/13_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('13_report_local.Rmd', output_file='$PWD/13_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else if (workflow.profile.contains('slurm'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else
        error "I'm not sure which to run, you must use local or slurm profiles"
    
}