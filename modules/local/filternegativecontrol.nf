process FILTERNEGATIVECONTROL{
    tag "$nc"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/decontam:1.2' : 
    'lorentzb/decontam:1.2' }"

    input:
    path('results')
    path('feature-table.tsv')
    path("metadata.tsv")
    path(control_script)
    val(nc)
 

    output:
    path("*.biom"), emit: filtered_table_biom
    path("*.tsv"), emit: filtered_table_tsv
    path("contam-features.tsv"), emit: contams_tsv
    

    script:

    '''
    #!/usr/bin/env bash

    echo ${nc} > nc_name.txt
    Rscript contam_script.r

    '''

}