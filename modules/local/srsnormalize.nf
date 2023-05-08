process SRSNORMALIZE{

    tag "$table"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/srs:1.0' : 
        'lorentzb/srs:1.0' }"

    input:

    path(table)
    val(srs_min)
    val(rare_val)
    
    output:
    path("*.tsv"), emit: tsv_normalized
    path("*.biom"), emit: biom_normalized
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript

    library(qiime2R)
    library(SRS)
    library(phyloseq)
    library(biomformat)

    #read in table from either decontam or results/qiime2/abundance_tables/feature-table.tsv
    un_rare_tab <- read.table('$table')


    if(file.exists('$srs_min')){
        srs_min <- as.numeric(readLines('$srs_min'))
    }


    if ($rare_val != 0){
        norm_tab <- SRS(un_rare_tab, $rare_val)
    } else {
        norm_tab <- SRS(un_rare_tab, srs_min)
    }

    rownames(norm_tab) <- rownames(un_rare_tab)

    # Save norm_tab.tsv

    #write.table(norm_tab, "normalized-table.tsv",row.names=F)
    write.table(norm_tab, "normalized-table.tsv")

    # convert norm_tab to biom 
    tmp <- make_biom(norm_tab)

    write_biom(tmp, "normalized-table.biom")

    """
}