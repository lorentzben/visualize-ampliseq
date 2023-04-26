// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process TSVTOQZA {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lorentzb/automate_16_nf:2.0':
        'lorentzb/automate_16_nf:2.0'}"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    
    tuple val(meta), path(tsv)
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
    
    biom add-metadata -i table.biom \
        -o md-table.biom \
        --observation-metadata-fp metadata.tsv

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
