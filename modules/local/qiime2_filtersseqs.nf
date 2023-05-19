process QIIME2_FILTERSEQS {
    tag "${filter}"
    label 'process_low'

    container "lorentzb/automate_16_nf:2.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(table)
    path(data) 
    path(taxonomy)
    val(filter)
    val(ioi)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: "--p-where \"[${ioi}]=\'${filter}\'\""
    def prefix = task.ext.prefix ?: "${filter}"
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"


    qiime feature-table filter-seqs \\
        --i-data ${data} \\
        --i-table ${table} \\
        --o-filtered-data ${prefix}_seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
