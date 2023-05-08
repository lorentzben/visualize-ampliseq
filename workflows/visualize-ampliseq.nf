#!/usr/bin/env nextflow


// Input 

if (params.input){
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    raw_tsv_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.tsv", checkIfExists: true)
    raw_biom_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.biom", checkIfExists: true)
    rooted_tree_ch = Channel.fromPath(params.input+"/qiime2/phylogenetic_tree/rooted-tree.qza", checkIfExists: true)
    
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

if(params.controls){
    controls_ch = Channel.fromPath(params.controls, checkIfExists:false)
    contam_script_ch = Channel.fromPath("${projectDir}/r_scripts/contam_script.r")
} else { controls_ch = Channel.empty() }

if(params.negative){
    nc_val_ch = Channel.of(params.negative)
} else { nc_val_ch = Channel.empty() }

if(params.mock){
    mock_val_ch = Channel.of(params.mock)
} else { mock_val_ch = Channel.empty() }

if(params.srs) {
    srs = true
    srs_curve_ch = Channel.fromPath("${projectDir}/r_scripts/srs_curve.rmd")
    srs_min_max_ch = Channel.fromPath("${projectDir}/python_scripts/my_count_table_min_max.py")
} else { 
    srs = false
    srs_curve_ch = Channel.empty()
    srs_min_max_ch = Channel.empty()
}

/*
/ Import Modules
*/
include { ORDERIOI } from "${projectDir}/modules/local/orderioi.nf"
include { CLEANUPRAWTSV; CLEANUPRAWTSV as CLEANUPFILTTSV; CLEANUPRAWTSV as CLEANUPFILTMOCKTSV; CLEANUPRAWTSV as CLEANUPFILTSRSTSV } from "${projectDir}/modules/local/cleanuprawtsv.nf"
include { CLEANUPRAWQZA } from "${projectDir}/modules/local/cleanuprawqza.nf"
include { TSVTOQZA; TSVTOQZA as TSVTOQZA2 } from "${projectDir}/modules/local/tsvtoqza.nf"
include { FILTERNEGATIVECONTROL } from "${projectDir}/modules/local/filternegativecontrol.nf"
include { QIIME2_FILTERSAMPLES as QIIME2_FILTERNC; QIIME2_FILTERSAMPLES as QIIME2_FILTERMOCK } from "${projectDir}/modules/local/qiime2_filtersamples.nf"
include { QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_NC; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_MOCK; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_CORE  } from "${projectDir}/modules/local/qiime2_export_absolute.nf"
include { SRSCURVE } from "${projectDir}/modules/local/srscurve.nf"
include { SRSNORMALIZE } from "${projectDir}/modules/local/srsnormalize.nf"


workflow VISUALIZEAMPLISEQ {
    //TODO see if this breaks it
    ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch
    ).ordered_ioi.set{ ord_ioi_ch }

    CLEANUPRAWTSV(raw_tsv_table_ch)

    CLEANUPRAWTSV.out.raw_table_tsv.set{ ch_raw_tsv_table }

    raw_mba_table = CLEANUPRAWTSV.out.raw_MbA_table_tsv

    

    CLEANUPRAWQZA(raw_biom_table_ch
    ).raw_table_qza.set { ch_raw_qza_table }

    ch_filtered_qza_table = Channel.empty()
    ch_filtered_tsv_table = Channel.empty()
    ch_normalized_qza = Channel.empty()
    ch_normalized_tsv = Channel.empty()

    ch_filtered_tsv_table.view()

    if(params.controls){
        filtered_table = FILTERNEGATIVECONTROL(input_ch, ch_raw_tsv_table, controls_ch, metadata_ch, contam_script_ch, nc_val_ch)
        tsv_map_1 = FILTERNEGATIVECONTROL.out.filtered_table_biom.map{
            it ->  [ [id: "Filtered-NC-Biom"], it ]
        }

        TSVTOQZA(tsv_map_1, metadata_ch
        ).qza.map{it.last()}.set{ ch_qza_filt_table }

        QIIME2_FILTERNC(metadata_ch, ch_qza_filt_table, nc_val_ch, ioi_ch
        ).qza.set{ ch_filtered_qza_table }

        QIIME2_EXPORT_ABSOLUTE_NC(QIIME2_FILTERNC.out.qza
        ).tsv.set { ch_messy_filtered_tsv_table }

        CLEANUPFILTTSV( ch_messy_filtered_tsv_table )

        CLEANUPFILTTSV.out.raw_table_tsv.set { ch_filtered_tsv_table }
        
        print("setting nc table to filtered!")
        ch_filtered_tsv_table.view()
        
        
    } 

    if(params.mock){
        mock_in_tsv = ch_filtered_tsv_table.ifEmpty(ch_raw_tsv_table)
        mock_in_qza = ch_filtered_qza_table.ifEmpty(ch_raw_qza_table)

        QIIME2_FILTERMOCK(metadata_ch, mock_in_qza, mock_val_ch, ioi_ch
            ).qza.set { ch_filtered_qza_table }
        QIIME2_EXPORT_ABSOLUTE_MOCK(QIIME2_FILTERMOCK.out.qza
            ).tsv.set { ch_messy_filtered_tsv_table }

        CLEANUPFILTMOCKTSV( ch_messy_filtered_tsv_table )

        CLEANUPFILTMOCKTSV.out.raw_table_tsv.set { ch_filtered_tsv_table }
        
        print("setting filtered to Mock!")
        ch_filtered_tsv_table.view()
        
    }

    if(srs){
        
        srs_in_tsv = ch_filtered_tsv_table.ifEmpty(ch_raw_tsv_table)
        srs_in_qza = ch_filtered_qza_table.ifEmpty(ch_raw_qza_table)

        SRSCURVE(srs_in_qza, srs_in_tsv, srs_curve_ch, srs_min_max_ch)

        SRSNORMALIZE(srs_in_tsv, SRSCURVE.out.min_val, params.rare
            ).tsv_normalized.set{ch_normalized_tsv}

        tsv_map_2 = SRSNORMALIZE.out.biom_normalized.map{
                it ->  [ [id: "SRS-Normalized-Biom"], it ]
        }

        TSVTOQZA2(tsv_map_2, metadata_ch
            ).qza.map{it.last()}.set{ch_normalized_qza}
        
        print("setting filtered to SRS")
        ch_filtered_tsv_table.view()
        
    } else{
        print("no normalization with SRS")
    }

    print('final table')
    ch_filtered_tsv_table.view()
}

    