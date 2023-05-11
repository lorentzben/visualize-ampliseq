process REPORT02GRAPHLANPHYLOGENETICTREE{

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/r_02:2.0' : 
        'lorentzb/r_02:2.0' }"

    tag "$reportName"
    label 'process_low'
   
    input:

    path("phylo_trees/*")
    path('item_of_interest.csv')
    path("02_report.Rmd")
    path("02_report_local.Rmd")

    output:

    path("02_report_*.html"), emit report02
     
    script:
    if(workflow.profile.contains('local'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('02_report_local.Rmd', output_file='$PWD/02_report_local_$dt.html', output_format='html_document',clean=TRUE,  knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('02_report_local.Rmd', output_file='$PWD/02_report_local_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else if (workflow.profile.contains('slurm'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.html', output_format='html_document',clean=TRUE,  knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else
        error "I'm not sure which to run, you must use local or slurm profiles"
    
}