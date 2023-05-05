#!/usr/bin/env nextflow

// Input 

if (params.input){
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
} else {
    log.error "Ampliseq input is required: please check params"
    System.exit(1)
}

if(params.metadata){
    metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
} else {
    log.error "Metadata is required: please check params"
    System.exit(1)
}

if(params.ioi){
    ioi_ch = Channel.of(params.ioi)
} else {
    log.error "Item of Interest is required: please check params"
    System.exit(1)
}

if(params.ordioi){
    ord_ioi_ch = Channel.fromPath(params.ordioi)
} else { ord_ioi_ch = Channel.empty() }

if(params.outdir){
    outdir_ch = Channel.fromPath(params.ordioi)
} else { outdir_ch = Channel.fromPath("results") }

if(params.rare){
    rare_val_ch = Channel.of(params.rare)
} else { rare_val_ch = Channel.of(0)}

if(params.control){
    controls_ch = Channel.fromPath(params.controls, checkIfExists:false)
} else { controls_ch = Channel.empty() }

if(params.negative){
    nc_val_ch = Channel.of(params.negative)
} else { nc_val_ch = Channel.empty() }

if(params.mock){
    mock_val_ch = Channel.of(params.mock)
} else { mock_val_ch = Channel.empty() }

if(params.srs) {
    srs = true
} else { srs = false }

/*
/ Import Modules
*/
include { ORDERIOI } from "${projectDir}/modules/local/orderioi.nf"
include { TSVTOQZA; TSVTOQZA as TSVTOQZA2 } from "${projectDir}/modules/local/tsvtoqza.nf"
include { QIIME2_FILTERSAMPLES as QIIME2_FILTERNC; QIIME2_FILTERSAMPLES as QIIME2_FILTERMOCK } from "${projectDir}/modules/local/qiime2_filtersamples.nf"
include { QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_NC; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_MOCK; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_CORE  } from "${projectDir}/modules/local/qiime2_export_absolute.nf"


//Not sure how to use this, but it's here
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

workflow VISUALIZEAMPLISEQ{
    //TODO see if this breaks it
    ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch
    ).ordered_ioi.set{ ord_ioi_ch }

}