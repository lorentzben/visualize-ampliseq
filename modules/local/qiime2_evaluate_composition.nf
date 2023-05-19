process QIIME2_EVALUATE_COMPOSITION {
    label 'process_low'

    container "lorentzb/automate_16_nf:2.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(experimentalTable)
    path(taxonomy)
    path(referenceTable)
    

    output:
    path("expected-observed-comparison.qzv"), emit: table_viz
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime taxa collapse \\
        --i-table ${experimentalTable} \\
        --i-taxonomy ${taxonomy} \\
        --p-level 7 \\
        --o-collapsed-table collapsed-table.qza

    qiime feature-table relative-frequency \\
        --i-table collapsed-table.qza \\
        --o-relative-frequency-table relative-table.qza

    # use quality control evaluate-seqs to check mock community
    qiime quality-control evaluate-composition \\
        --i-expected-features ${referenceTable} \\
        --i-observed-features relative-table.qza \\
        --o-visualization expected-observed-comparison.qzv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
