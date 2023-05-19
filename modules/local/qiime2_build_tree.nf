process QIIME2_BUILD_ROOTED_TREE {
    label 'process_low'

    container "lorentzb/automate_16_nf:2.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(repSeqs)
     

    output:
    path("aligned-rep-seqs.qza"), emit: alignment
    path("masked-aligned-rep-seqs.qza"), emit: maskedAlignment
    path("unrooted-tree.qza"), emit: unrootedTree
    path("rooted-tree.qza"), emit: rootedTree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    def prefix = task.ext.prefix ?: "${repSeqs}"
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences ${repSeqs} \\
        --o-alignment aligned-rep-seqs.qza \\
        --o-masked-alignment masked-aligned-rep-seqs.qza \\
        --o-tree unrooted-tree.qza \\
        --o-rooted-tree rooted-tree.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
