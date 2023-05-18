process QIIME2_FILTERSAMPLES {
    tag "${filter}"
    label 'process_low'

    container "lorentzb/automate_16_nf:2.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(metadata)
    path(table) 
    val(filter)
    val(ioi)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "--p-where \"[${ioi}]!=\'${filter}\'\""
    def prefix = task.ext.prefix ?: "${filter}"
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime feature-table filter-samples \\
        --i-table ${table} \\
        --m-metadata-file ${metadata} \\
        $args \\
        --o-filtered-table ${prefix}.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
