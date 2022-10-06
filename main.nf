#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.analysis = "${projectDir}/ampliseq_results"
params.outdir = "results"

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         analysis : ${params.analysis}
         outdir   : ${params.outdir}
         """
         .stripIndent()


process REPORT01BARPLOT{
    input:

    path analysis_files

    script:

    """
    #!/usr/bin/env python3

    print(${analysis_files})    
    """

}

workflow{
    input_ch = Channel.fromPath(params.analysis)
    REPORT01BARPLOT(input_ch)
}