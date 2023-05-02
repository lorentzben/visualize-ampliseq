process QIIME2_EXPORT_ABSOLUTE {
    label 'process_low'

    container "quay.io/qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(table)

    output:
    path("feature-table.tsv")        , emit: tsv
    path("feature-table.biom")       , emit: biom
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export \\
        --input-path ${table} \\
        --output-path table
    cp table/feature-table.biom .

    #produce raw count table "table/feature-table.tsv"
    biom convert \\
        -i table/feature-table.biom \\
        -o feature-table.tsv \\
        --to-tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
