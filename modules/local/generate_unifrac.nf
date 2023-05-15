process GENERATEUNIFRAC{

    tag "$distances"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/automate_16_nf:2.0' : 
        'lorentzb/automate_16_nf:2.0' }"

    input:

    path(distances)
    path("metadata.tsv")
    file('item_of_interest.csv')

    output:

    path("*-pairwise.tsv"), emit: pairwise
     
    script:

    '''
    #!/usr/bin/env bash
    IOI=$(cat item_of_interest.csv)

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    qiime diversity beta-group-significance --i-distance-matrix weighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column $IOI \
    --p-pairwise \
    --o-visualization weighted-unifrac.qzv

    qiime tools export \
    --input-path weighted-unifrac.qzv \
    --output-path weighted-unifrac

    cp weighted-unifrac/raw_data.tsv ./weighted-unifrac-pairwise.tsv

    qiime diversity beta-group-significance --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column $IOI \
    --p-pairwise \
    --o-visualization unweighted-unifrac.qzv

    qiime tools export \
    --input-path unweighted-unifrac.qzv \
    --output-path unweighted-unifrac

    cp  unweighted-unifrac/raw_data.tsv ./unweighted-unifrac-pairwise.tsv
    '''
}