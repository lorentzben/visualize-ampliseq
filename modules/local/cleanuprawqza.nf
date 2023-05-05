process CLEANUPRAWQZA{

    tag "$table_biom"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/automate_16_nf:2.0' : 
    'lorentzb/automate_16_nf:2.0' }"

    input:

    path (table_biom)
    //results/qiime2/abundance_tables/feature-table.biom

    output:
    
    path ("raw-feature-table.qza"), emit: raw_table_qza

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time
    import os

    create_qza_command = "qiime tools import \
    --input-path $table_biom  \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path raw-feature-table.qza"
    result = subprocess.run([create_qza_command], shell=True)
    """
}