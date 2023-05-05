#!/usr/bin/env nextflow

// Input 

if (params.input){
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
} else {
    log.error "Ampliseq input is required: please check params"
    System.exit(1)
}




params.metadata = "${projectDir}/metadata.tsv"
params.ioi = "treatment"
params.ordioi = "ordered_item_of_interest.csv"
params.outdir = "results"
params.rare = 0
params.controls = ""
params.srs = false
params.mock = ""
params.negative = ""

//Not sure how to use this, but it's here
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

workflow VISUALIZEAMPLISEQ{


}