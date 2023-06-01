process GENERATEBIOMFORGRAPHLAN{
    
    tag "$ioi"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/automate_16_nf:2.0' : 
    'lorentzb/automate_16_nf:2.0' }"

    input:

    path(metadata)
    val(ioi)
    path("filter_samples.py")
    path("taxonomy.qza")
    path("feature-table.qza")
    path("feature-table.tsv")
    val(mock)
    val(nc)
    

    output:
    
    path("biom_tabs/*"), emit: graphlan_biom
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time
    import os

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    
    ioi = '${ioi}'

    # TODO check that all ioi set present in table?

    feature_table = pd.read_table('feature-table.tsv', sep=' ')
    samples = list(feature_table.columns)

    metadata_table = metadata_table[metadata_table['ID'].isin(samples)]

    ioi_set = set(metadata_table[\"${ioi}\"])
    
    ioi_set.discard("$mock")
    ioi_set.discard("$nc")
        


    subprocess.run(['mkdir phylo_trees'], shell=True)
    subprocess.run(['mkdir biom_tabs'], shell=True)

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi
        
        filter_command = "python3 filter_samples.py -m ${metadata} -i ${ioi} -c "+str(item)
        result = subprocess.run([filter_command], shell=True)

        time.sleep(2)

        # adds taxonomic info needed for plotting
        collapse_command = 'qiime taxa collapse \
        --i-table '+str(item)+'-filtered-table.qza \
        --o-collapsed-table collapse-'+str(item)+'-table.qza \
        --p-level 7 \
        --i-taxonomy taxonomy.qza'

        result = subprocess.run([collapse_command], shell=True)

        # exports artifact so that the next step can collect it
        export_command='qiime tools export \
        --input-path collapse-'+str(item)+'-table.qza \
        --output-path collapse-'+str(item)+'-frequency/'
        
        result = subprocess.run([export_command], shell=True)

        # turns feature table into a human-reable format
        biom_command = 'biom convert -i collapse-'+str(item)+\
        '-frequency/feature-table.biom -o otu-'+str(item)+\
        '-table.tsv --to-tsv --header-key taxonomy'

        result = subprocess.run([biom_command], shell=True)

        # formatting the table so that it is in the correct order
        table = pd.read_table(\"otu-"+str(item)+"-table.tsv\", sep='\t', header=1)
        table = table.drop(columns=['taxonomy'])
        table = table.rename(columns={'#OTU ID':'taxonomy'})
        tax = table.pop('taxonomy')
        insertion_site = len(table.columns)
        table.insert(insertion_site, 'taxonomy', tax)
        table.insert(0, 'OTU_ID', np.arange(len(table)))
        table.to_csv('otu-'+str(item)+'-mod-table.tsv', sep='\t', index=False)

        # human readable table into compressed computer-readble format
        biom_format_command='biom convert -i otu-'+str(item)+ \
        '-mod-table.tsv -o '+str(item)+'-otu-table-mod.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy'

        result = subprocess.run([biom_format_command], shell=True)

        result = subprocess.run(['cp '+str(item)+'-otu-table-mod.biom biom_tabs'],shell=True)
    """

}