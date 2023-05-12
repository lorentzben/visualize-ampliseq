process GENERATERAREFACTIONCURVE{

    tag "$rare_val"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/automate_16_nf:2.0' : 
        'lorentzb/automate_16_nf:2.0' }"

    input:
    
    path('metadata.tsv')
    path('table.qza')
    path('rooted-tree.qza')
    val(rare_val)
    path("table.tsv")

    output:
    path("rarefact/observed_features.csv"), emit: rareVector

    script:
    def table_tsv = table_tsv.name != 'NO_FILE' ? "$table_tsv" : ''

    """
    #!/usr/bin/env python3


    from qiime2.plugins.diversity.visualizers import alpha_rarefaction
    from qiime2.plugins import diversity
    from qiime2 import Metadata
    from qiime2.plugins import feature_table
    from qiime2 import Artifact
    import pandas as pd
    import sys
    import warnings
    import os

    table = Artifact.load('table.qza')
    metadata = Metadata.load('metadata.tsv')
    rooted_tree = Artifact.load('rooted-tree.qza')

    # if the default value aka use count_table_minmax_reads
    if $rare_val == 0:
        uncompress_table='table.tsv'

        # adapted from count_table_minmax_reads.py @author Daniel Straub
        # collected from nf-core/ampliseq
        # read tsv and skip first two rows
        data = pd.read_csv('table.tsv', sep="\t", skiprows=[0], header=None)  # count table

        # drop feature ids
        df = data.drop(data.columns[0], axis=1)

        # make sums
        sums = df.sum()

        # we want minimum values
        mindepth = int(sums.min())
        maxdepth = int(sums.max())

        if mindepth > 10000:
            print("Use the sampling depth of " +str(mindepth)+" for rarefaction")
        elif mindepth < 10000 and maxdepth > 5000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is quite small for rarefaction")
        elif mindepth < 5000 and mindepth > 1000:
            print("WARNING The sampling depth of "+str(mindepth)+" is very small for rarefaction")
        elif mindepth < 1000:
            print("WARNING The sampling depth of "+str(mindepth)+" seems too small for rarefaction")
        else:
            print("ERROR this shouldn't happen")
            exit(1)

        # TODO convert this from bash to python
        #check values
        
        if maxdepth > 75000:
            maxdepth = 75000
        
        if maxdepth > 5000:
            maxsteps=250
        else:
            maxsteps=(maxdepth/20)

        rarefact = alpha_rarefaction(table=table, max_depth=maxdepth, phylogeny=rooted_tree, steps=maxsteps)
        file = open("rarefaction.txt", "w")
        file.write(str(mindepth))
        file.close 
    
    # else if user submits the rarefaction depth they want to use based on rarefaction plot
    else: 
        rarefact = alpha_rarefaction(table=table, max_depth=$rare_val, phylogeny=rooted_tree)
        file = open("rarefaction.txt", "w")
        file.write(str($rare_val))
        file.close 

    
    Artifact.export_data(rarefact.visualization,'rarefact')


    """
}