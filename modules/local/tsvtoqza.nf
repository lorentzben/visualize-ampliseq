process TSVTOQZA {
    tag "$meta.id"
    label 'process_low'

    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lorentzb/automate_16_nf:2.0':
        'lorentzb/automate_16_nf:2.0'}"

    input:

    tuple val(meta), path(biom)
    path "metadata.tsv"

    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env bash
    
    biom add-metadata -i $biom -o md-table.biom --observation-metadata-fp metadata.tsv

    qiime tools import \
        --input-path md-table.biom \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path feature-table.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biom: \$(echo \$(biom --version 2>&1))
        qiime: \$(echo \$(qiime info 2>&1))
    END_VERSIONS
    """
}
