process LEFSEFORMAT{

    tag "$ioi"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/qiime2lefse:1.0' : 
        'lorentzb/qiime2lefse:1.0' }"

    input:
    val(ioi)
    path("feature-table.qza")
    path("rooted-tree.qza")
    path("taxonomy.qza")
    path("metadata.tsv")
    path("qiime_to_lefse.R")
      
    output:
    path ("combos/*"), emit : combos

    script:
    """
    #!/usr/bin/env bash
    mkdir combos
    Rscript qiime_to_lefse.R ${ioi}
    mv lefse_formatted.txt combos/
    """
}


process LEFSEANALYSIS{

    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/py2_env:1.0' : 
        'lorentzb/py2_env:1.0' }"

    input:

    path("combos/*")
    path("lefse_analysis.sh")
    path(plot_clado)
    path(plot_res)

    output:
    path("lefse_images/*"), emit: lefse_images

    script:
    """
    #!/usr/bin/env bash
    mkdir lefse_images
    bash lefse_analysis.sh
    """
}