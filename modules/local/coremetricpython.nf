process COREMETRICPYTHON{

    tag "$rare_val"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/automate_16_nf:2.0' : 
    'lorentzb/automate_16_nf:2.0' }"

    input:

    path(metadata)
    path(table_qza)
    path(table_tsv)
    path('rooted-tree.qza')
    val(rare_val)
    //TODO if we want to use a raw table for rare_val=0 
    //path(raw_table)
    

    output:

    path("diversity_core/*_pcoa_results.qza")   , emit: pcoa
    path("diversity_core/*_vector.qza")         , emit: vector
    path("diversity_core/*_distance_matrix.qza"), emit: distance
    path("*rarefaction.txt") , emit: depth
    path("diversity_core/rarefied_table.qza"), emit: rare_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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

  
    unrarefied_table = Artifact.load('$table_qza')

    os.mkdir("diversity_core")

    warnings.filterwarnings('ignore')

    metadata = Metadata.load('$metadata')

    rooted_tree = Artifact.load('rooted-tree.qza')

    # if the default value aka use count_table_minmax_reads
    if $rare_val == 0:

        # adapted from count_table_minmax_reads.py @author Daniel Straub
        # collected from nf-core/ampliseq
        # read tsv and skip first two rows
        if "$table_tsv".find("normalized") >= 0 :
            data = pd.read_csv('$table_tsv', sep="\s", skiprows=[0], header=None)  # count table
        else:
            data = pd.read_csv('$table_tsv', sep="\t", skiprows=[0], header=None)  # count table

        # drop feature ids
        df = data.drop(data.columns[0], axis=1)

        # make sums
        sums = df.sum()

        # we want minimum values
        mindepth = int(sums.min())

        if mindepth > 10000:
            print("Use the sampling depth of " +str(mindepth)+" for rarefaction")
        elif mindepth < 10000 and mindepth > 5000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is quite small for rarefaction")
        elif mindepth < 5000 and mindepth > 1000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is very small for rarefaction")
        elif mindepth < 1000: 
            print("WARNING The sampling depth of "+str(mindepth)+" seems too small for rarefaction")
        else:
            print("ERROR this shouldn't happen")
            exit(1)

        core = diversity.pipelines.core_metrics_phylogenetic(unrarefied_table, rooted_tree, mindepth, metadata)
        file = open("rarefaction.txt", "w")
        file.write(str(mindepth))
        file.close 
    
    # else if user submits the rarefaction depth they want to use based on rarefaction plot
    else: 
        core = diversity.pipelines.core_metrics_phylogenetic(unrarefied_table, rooted_tree, $rare_val, metadata)
        file = open("rarefaction.txt", "w")
        file.write(str($rare_val))
        file.close 

    Artifact.save(core[0], "diversity_core/rarefied_table")
    Artifact.save(core[1], "diversity_core/faith_pd_vector")
    Artifact.save(core[2], "diversity_core/observed_features_vector")
    Artifact.save(core[3], "diversity_core/shannon_vector")
    Artifact.save(core[4], "diversity_core/evenness_vector")
    Artifact.save(core[5], "diversity_core/unweighted_unifrac_distance_matrix")
    Artifact.save(core[6], "diversity_core/weighted_unifrac_distance_matrix")
    Artifact.save(core[7], "diversity_core/jaccard_distance_matrix")
    Artifact.save(core[8], "diversity_core/bray_curtis_distance_matrix")
    Artifact.save(core[9], "diversity_core/unweighted_unifrac_pcoa_results")
    Artifact.save(core[10], "diversity_core/weighted_unifrac_pcoa_results")
    Artifact.save(core[11], "diversity_core/jaccard_pcoa_results")
    Artifact.save(core[12], "diversity_core/bray_curtis_pcoa_results")
    Artifact.save(core[13], "diversity_core/unweighted_unifrac_emperor")
    Artifact.save(core[14], "diversity_core/weighted_unifrac_emperor")
    Artifact.save(core[15], "diversity_core/jaccard_emperor")
    Artifact.save(core[16], "diversity_core/bray_curtis_emperor")  

    """
}