process RUNGRAPHLAN{
   
    tag "$ioi"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/py2_test:2.0' : 
        'lorentzb/py2_test:2.0' }"

    input:

    path(metadata)
    val(ioi)
    path("taxonomy.qza")
    path("graph.sh")
    path("biom_tabs/*")
   
    output:
    path("phylo_trees/*"), emit: graphlan_dir

    script:
    """
    #!/usr/bin/env python2
    import subprocess
    import csv
    import pandas as pd
    import numpy as np 
    import time
    import os

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    ioi_set = set(metadata_table[\"${ioi}\"])
    ioi = '${ioi}'

    os.system('cp biom_tabs/*-otu-table-mod.biom .')

    os.system('mkdir phylo_trees')

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi

        # Outputs the current ioi so that it can be annotated in the graphlan image
        with open('current.txt', 'w') as file:
            file.write(item)

        # bash script call to handle the steps within a conda python 2.7.17 envionment
        generate_image_command = 'bash graph.sh'
        result = os.system(generate_image_command)

    rename_image = 'cp *_image_graph.png phylo_trees/.'
    result = os.system(rename_image)

    rename_pdf_image = 'cp *_image_pdf_graph.png phylo_trees/.'
    result = os.system(rename_pdf_image)


    """

}