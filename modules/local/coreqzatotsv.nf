process COREQZATOTSV{

    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/automate_16_nf:2.0' : 
        'lorentzb/automate_16_nf:2.0' }"

    input:

    path (diversity)
    
    output:

    path("*_vector.tsv"), emit: vector

    script:

    """
    #!/usr/bin/env python3

    from qiime2.plugins import diversity
    from qiime2 import Metadata
    from qiime2.plugins import feature_table
    from qiime2 import Artifact
    import pandas as pd
    import sys
    import warnings
    import os

    diversity_names = '$diversity'
    diversity_names = diversity_names.split(' ')

    for item in diversity_names:

        diversity_obj = Artifact.load(item)

        Artifact.export_data(diversity_obj,'.')
        artifact_name = item
        filename = str(artifact_name.split('.')[0]+'.tsv')

        os.rename('alpha-diversity.tsv', filename)
    """
}