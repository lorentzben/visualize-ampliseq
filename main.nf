#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { VISUALIZEAMPLISEQ } from "${projectDir}/workflows/visualize-ampliseq.nf"

workflow VISUALIZE_AMPLISEQ{
    VISUALIZEAMPLISEQ()
}

workflow {
    VISUALIZE_AMPLISEQ()
}