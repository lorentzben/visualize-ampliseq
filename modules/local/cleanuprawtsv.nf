process CLEANUPRAWTSV{
    tag "$table"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/srs:1.0' : 
        'lorentzb/srs:1.0' }"

    input:
    
    path(table)
 

    output:
    path('raw_table.tsv'), emit: raw_table_tsv
    path('raw_table_MbA.tsv'), emit: raw_MbA_table_tsv
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    
    new_table <- read.csv('$table',sep='\t', skip=1)
    new_table_df <- data.frame(new_table)
    colnames(new_table_df)[1] <- "#NAME"
    rownames(new_table_df) <- new_table_df[,1]
    write.table(new_table_df, "raw_table_MbA.tsv",row.names=F,sep='\t')
    collen <- dim(new_table_df)[2]
    new_table_df <- new_table_df[,2:collen]
    write.table(new_table_df, "raw_table.tsv", row.names=T ,sep='\t')
    system("rm $table")
    """

}