#!/usr/bin/env nextflow


// Input 

if (params.input){
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    raw_tsv_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.tsv", checkIfExists: true)
    raw_biom_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.biom", checkIfExists: true)
    rooted_tree_ch = Channel.fromPath(params.input+"/qiime2/phylogenetic_tree/rooted-tree.qza", checkIfExists: true)
    asv_tsv_ch = Channel.fromPath(params.input+"/dada2/ASV_tax_species.tsv", checkIfExists: true)
    ch_overall_summary = Channel.fromPath(params.input+"/overall_summary.tsv", checkIfExists: true)
    ch_tree_nwk = Channel.fromPath(params.input+"/qiime2/phylogenetic_tree/tree.nwk")
    filter_samples_ch = Channel.fromPath("${projectDir}/python_scripts/filter_samples.py")
    
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

if(params.refSeq){
    refrence_seq_ch = Channel.fromPath(params.refSeq)
    rep_seq_ch = Channel.fromPath(params.input+"/qiime2/representative_sequences/filtered-sequences.qza", checkIfExists: true)
} else { 
    refrence_seq_ch = Channel.empty() 
    rep_seq_ch = Channel.empty()
}

if(params.refTab){
    reference_table_ch = Channel.fromPath(params.refTab)
} else{ reference_table_ch = Channel.empty() }


if(params.srs) {
    srs = true
    srs_curve_ch = Channel.fromPath("${projectDir}/r_scripts/srs_curve.rmd")
    srs_min_max_ch = Channel.fromPath("${projectDir}/python_scripts/my_count_table_min_max.py")
} else { 
    srs = false
    srs_curve_ch = Channel.empty()
    srs_min_max_ch = Channel.empty()
    
}

if(params.report){
    graph_sh_ch = Channel.fromPath("${projectDir}/bash_scripts/graph.sh")
    qiime_to_lefse_ch = Channel.fromPath("${projectDir}/r_scripts/qiime_to_lefse.R")
    lefse_analysis_ch = Channel.fromPath("${projectDir}/bin/lefse_analysis.sh")
    plot_clado_file_ch = Channel.fromPath("${projectDir}/python_scripts/plot_cladogram.py")
    plot_res_file_ch = Channel.fromPath("${projectDir}/python_scripts/plot_res.py")

    report_one_ch = Channel.fromPath("${projectDir}/report_gen_files/01_report_MbA.Rmd")
    report_two_ch = Channel.fromPath("${projectDir}/report_gen_files/02_report.Rmd")
    report_two_local_ch = Channel.fromPath("${projectDir}/report_gen_files/02_report_local.Rmd")
    report_three_ch = Channel.fromPath("${projectDir}/report_gen_files/03_report.Rmd")
    report_four_ch = Channel.fromPath("${projectDir}/report_gen_files/04_report.Rmd")
    report_five_ch = Channel.fromPath("${projectDir}/report_gen_files/05_report.Rmd")
    report_six_ch = Channel.fromPath("${projectDir}/report_gen_files/06_report.Rmd")
    report_six_b_ch = Channel.fromPath("${projectDir}/report_gen_files/06b_report.Rmd")
    report_seven_ch = Channel.fromPath("${projectDir}/report_gen_files/07_report.Rmd")
    report_eight_ch = Channel.fromPath("${projectDir}/report_gen_files/08_report.Rmd")
    report_nine_ch = Channel.fromPath("${projectDir}/report_gen_files/09_report.Rmd")
    report_ten_ch = Channel.fromPath("${projectDir}/report_gen_files/10_report.Rmd")
    report_eleven_ch = Channel.fromPath("${projectDir}/report_gen_files/11_report.Rmd")
    report_twelve_ch = Channel.fromPath("${projectDir}/report_gen_files/12_report.Rmd")
    report_thirteen_ch = Channel.fromPath("${projectDir}/report_gen_files/13_report.Rmd")
    report_thirteen_local_ch = Channel.fromPath("${projectDir}/report_gen_files/13_report_local.Rmd")
    report_fourteen_ch = Channel.fromPath("${projectDir}/report_gen_files/14_report.Rmd")
}

/*
/ Import Modules
*/
include { ORDERIOI } from "${projectDir}/modules/local/orderioi.nf"
include { REFORMATANDQZATAX } from "${projectDir}/modules/local/reformatandqzatax.nf"
include { CLEANUPRAWTSV; CLEANUPRAWTSV as CLEANUPFILTTSV; CLEANUPRAWTSV as CLEANUPFILTMOCKTSV; CLEANUPRAWTSV as CLEANUPFILTSRSTSV; CLEANUPRAWTSV as CLEANUPNORMTSV } from "${projectDir}/modules/local/cleanuprawtsv.nf"
include { CLEANUPRAWQZA } from "${projectDir}/modules/local/cleanuprawqza.nf"
include { TSVTOQZA; TSVTOQZA as TSVTOQZA2 } from "${projectDir}/modules/local/tsvtoqza.nf"
include { FILTERNEGATIVECONTROL } from "${projectDir}/modules/local/filternegativecontrol.nf"
include { QIIME2_FILTERSAMPLES as QIIME2_FILTERNC; QIIME2_FILTERSAMPLES as QIIME2_FILTERMOCK } from "${projectDir}/modules/local/qiime2_filternotsamples.nf"
include { QIIME2_FILTERONLYSAMPLES as QIIME2_ONLYMOCK } from "${projectDir}/modules/local/qiime2_filtersamples.nf"
include { QIIME2_FILTERSEQS; QIIME2_FILTERSEQS as QIIME2_FILTER_REPSEQS } from "${projectDir}/modules/local/qiime2_filtersseqs.nf"
include { QIIME2_BUILD_ROOTED_TREE } from "${projectDir}/modules/local/qiime2_build_tree.nf"
include { QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_NC; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_MOCK; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_CORE  } from "${projectDir}/modules/local/qiime2_export_absolute.nf"
include { SRSCURVE } from "${projectDir}/modules/local/srscurve.nf"
include { SRSNORMALIZE } from "${projectDir}/modules/local/srsnormalize.nf"
include { GENERATEBIOMFORGRAPHLAN } from "${projectDir}/modules/local/generatebiomforgraphlan.nf"
include { RUNGRAPHLAN } from "${projectDir}/modules/local/graphlan.nf"
include { COREMETRICPYTHON } from "${projectDir}/modules/local/coremetricpython.nf"
include { COREQZATOTSV } from "${projectDir}/modules/local/coreqzatotsv.nf"
include { GENERATERAREFACTIONCURVE } from "${projectDir}/modules/local/generate_rarefaction_curve.nf"
include { GENERATEUNIFRAC } from "${projectDir}/modules/local/generate_unifrac.nf"
include { LEFSEFORMAT; LEFSEANALYSIS } from "${projectDir}/modules/local/lefse.nf"
include { QIIME2_EVALUATE_SEQS } from "${projectDir}/modules/local/qiime2_evaluate_seqs.nf"
include { QIIME2_EVALUATE_COMPOSITION } from "${projectDir}/modules/local/qiime2_evaluate_composition.nf"
include { REPORT01BARPLOT } from "${projectDir}/modules/local/renderreport01.nf"
include { REPORT02GRAPHLANPHYLOGENETICTREE } from "${projectDir}/modules/local/renderreport02.nf"
include { REPORT03HEATMAP } from "${projectDir}/modules/local/renderreport03.nf"
include { REPORT04ALPHATABLE } from "${projectDir}/modules/local/renderreport04.nf"
include { REPORT05ALPHABOXPLOT } from "${projectDir}/modules/local/renderreport05.nf"
include { REPORT06ORDINATION; REPORT06BNMDSORDINATION } from "${projectDir}/modules/local/renderreport06.nf"
include { REPORT07RAREFACTION } from "${projectDir}/modules/local/renderreport07.nf"
include { REPORT08RANKEDABUNDANCE } from "${projectDir}/modules/local/renderreport08.nf"
include { REPORT09UNIFRACHEATMAP } from "${projectDir}/modules/local/renderreport09.nf"
include { REPORT10BETABOXPLOT } from "${projectDir}/modules/local/renderreport10.nf"
include { REPORT11UPGMA } from "${projectDir}/modules/local/renderreport11.nf"
include { REPORT12PERMANOVA } from "${projectDir}/modules/local/renderreport12.nf"
include { REPORT13LEFSE } from "${projectDir}/modules/local/renderreport13.nf"
include { REPORT14CITATIONS } from "${projectDir}/modules/local/renderreport14.nf"

workflow VISUALIZEAMPLISEQ {
    ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch
    ).ordered_ioi.set{ ord_ioi_ch }

    CLEANUPRAWTSV(raw_tsv_table_ch)

    CLEANUPRAWTSV.out.raw_table_tsv.set{ ch_raw_tsv_table }

    raw_mba_table = CLEANUPRAWTSV.out.raw_MbA_table_tsv

    REFORMATANDQZATAX(asv_tsv_ch
        ).tax_qza.set{ ch_tax_qza }

    CLEANUPRAWQZA(raw_biom_table_ch
    ).raw_table_qza.set { ch_raw_qza_table }

    ch_filtered_qza_table = Channel.empty()
    ch_filtered_tsv_table = Channel.empty()
    ch_normalized_qza = Channel.empty()
    ch_normalized_tsv = Channel.empty()


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
        
    } 

    if(params.mock){
        if(params.controls){
            // Yes Mock Yes Negative Control
            mock_in_tsv = ch_filtered_tsv_table
            mock_in_qza = ch_filtered_qza_table

        } else{
            // Yes Mock No Negative Control
            mock_in_tsv = ch_raw_tsv_table
            mock_in_qza = ch_raw_qza_table
        }
        

        QIIME2_FILTERMOCK(metadata_ch, mock_in_qza, mock_val_ch, ioi_ch
            ).qza.set { ch_filtered_qza_table }
        QIIME2_EXPORT_ABSOLUTE_MOCK(QIIME2_FILTERMOCK.out.qza
            ).tsv.set { ch_messy_filtered_tsv_table }

        CLEANUPFILTMOCKTSV( ch_messy_filtered_tsv_table )

        CLEANUPFILTMOCKTSV.out.raw_table_tsv.set { ch_filtered_tsv_table }
        
    }

    if(srs){
        if(params.controls){
            if(params.mock){
                // Yes Negative Control Yes Mock Yes SRS
                srs_in_tsv = ch_filtered_tsv_table
                srs_in_qza = ch_filtered_qza_table

            } else{
                // Yes Negative Control No Mock Yes SRS
                srs_in_tsv = ch_filtered_tsv_table
                srs_in_qza = ch_filtered_qza_table
            }
        } else if(params.mock){
            // No Negative Control Yes Mock Yes SRS
            srs_in_tsv = ch_filtered_tsv_table
            srs_in_qza = ch_filtered_qza_table

        } else{
            // No Negative Control No Mock Yes SRS
            srs_in_tsv = ch_raw_tsv_table
            srs_in_qza = ch_raw_qza_table

        }
        
        SRSCURVE(srs_in_qza, srs_in_tsv, srs_curve_ch, srs_min_max_ch)

        SRSNORMALIZE(srs_in_tsv, SRSCURVE.out.min_val, params.rare
            ).tsv_normalized.set{ch_normalized_tsv}

        tsv_map_2 = SRSNORMALIZE.out.biom_normalized.map{
                it ->  [ [id: "SRS-Normalized-Biom"], it ]
        }

        TSVTOQZA2(tsv_map_2, metadata_ch
            ).qza.map{it.last()}.set{ch_normalized_qza}
        
        
        
    } else{
        //print("no normalization with SRS")
    }

    if (params.controls){
        if(params.mock){
            if(srs){
                // Yes NC Yes Mock Yes SRS
                final_table_tsv = ch_normalized_tsv
                final_table_qza = ch_normalized_qza

            } else {
                // Yes NC Yes Mock No SRS
                final_table_tsv = ch_filtered_tsv_table
                final_table_qza = ch_filtered_qza_table
            }
            
        } else {
            if(srs){
                // Yes NC No Mock Yes SRS
                final_table_tsv = ch_normalized_tsv
                final_table_qza = ch_normalized_qza
            } else {
                // Yes NC No Mock No SRS
                final_table_tsv = ch_filtered_tsv_table
                final_table_qza = ch_filtered_qza_table
            }
            
        }
    } else {
        if(params.mock){
            if(srs){
                // No NC Yes Mock Yes SRS
                final_table_tsv = ch_normalized_tsv
                final_table_qza = ch_normalized_qza
            } else {
                // No NC Yes Mock No SRS
                final_table_tsv = ch_filtered_tsv_table
                final_table_qza = ch_filtered_qza_table
            }
            
        } else {
            if(srs){
                // No NC No Mock Yes SRS
                final_table_tsv = ch_normalized_tsv
                final_table_qza = ch_normalized_qza
            } else {
                // No NC No Mock No SRS
                final_table_tsv = ch_raw_tsv_table
                final_table_qza = ch_raw_qza_table
            }
        }
    }

    GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, filter_samples_ch, ch_tax_qza, final_table_qza, nc_val_ch.ifEmpty("N/A"), mock_val_ch.ifEmpty("N/A")
        ).graphlan_biom.set{ ch_graphlan_biom }

    RUNGRAPHLAN(metadata_ch, ioi_ch, ch_tax_qza, graph_sh_ch, ch_graphlan_biom
        ).graphlan_dir.set{ ch_graphlan_dir }
    
    //TODO Use the filtered table to filter repseqs

    QIIME2_FILTER_REPSEQS(final_table_qza, rep_seq_ch
        ).qza.set{ ch_new_rep_seq }

    //TODO Rebuild the rooted_tree

    QIIME2_BUILD_ROOTED_TREE( rooted_tree_ch, ch_new_rep_seq 
        ).rootedTree{ new_rooted_tree_ch }

    COREMETRICPYTHON(metadata_ch, final_table_qza, final_table_tsv, new_rooted_tree_ch, rare_val_ch
        ).rare_table.set{ ch_norm_qza_table }

    COREMETRICPYTHON.out.pcoa.set{ ch_core_pcoa }
    COREMETRICPYTHON.out.vector.set{ ch_core_vector }
    COREMETRICPYTHON.out.distance.set{ ch_core_distance }

    QIIME2_EXPORT_ABSOLUTE_CORE(ch_norm_qza_table
        ).tsv.set{ ch_norm_messy_tsv_table }
    
    CLEANUPNORMTSV( ch_norm_messy_tsv_table 
        ).raw_table_tsv.set{ ch_norm_tsv_table }

    CLEANUPNORMTSV.out.raw_MbA_table_tsv.set{ ch_norm_MBA_tsv_table }
    COREQZATOTSV(COREMETRICPYTHON.out.vector
        ).vector.set{ ch_core_vector_tsv }
    
    GENERATEUNIFRAC(ch_core_distance, metadata_ch, ioi_ch
        ).pairwise.set{ unifrac_pairwise_ch }

    LEFSEFORMAT(ioi_ch, ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, qiime_to_lefse_ch
        ).combos.set{ ch_lefse_combos }
    LEFSEANALYSIS( ch_lefse_combos, lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch
        ).lefse_images.set{ ch_lefse_images }

    REPORT01BARPLOT("Report_01", asv_tsv_ch, ch_tree_nwk, metadata_ch, report_one_ch, ioi_ch, ch_norm_MBA_tsv_table, ch_norm_qza_table)
    REPORT02GRAPHLANPHYLOGENETICTREE( "Report_02", ch_graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
    REPORT03HEATMAP("Report_03", ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi_ch, ch_overall_summary)
    REPORT04ALPHATABLE("Report_04", ch_core_vector_tsv, ioi_ch, report_four_ch)
    REPORT05ALPHABOXPLOT("Report_05", ch_core_vector_tsv, ioi_ch, ord_ioi_ch, metadata_ch, report_five_ch)
    REPORT06ORDINATION("Report_06", ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, ioi_ch, ord_ioi_ch, report_six_ch, ch_core_pcoa, ch_core_vector)
    REPORT06BNMDSORDINATION("Report_06b",ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, ioi_ch, ord_ioi_ch, report_six_b_ch, ch_core_pcoa, ch_core_vector)

    if (!params.srs){
        // Note will always take a raw table
        GENERATERAREFACTIONCURVE(metadata_ch, ch_raw_qza_table, new_rooted_tree_ch, rare_val_ch, ch_raw_tsv_table
            ).rareVector.set{ ch_rare_vector }
        REPORT07RAREFACTION("Report_07", ioi_ch, ord_ioi_ch, report_seven_ch, ch_rare_vector, metadata_ch)
    }

    REPORT08RANKEDABUNDANCE("Report_08", ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, ioi_ch, ord_ioi_ch, report_eight_ch)
    REPORT09UNIFRACHEATMAP("Report_09", ioi_ch, ord_ioi_ch, metadata_ch, ch_core_distance, report_nine_ch)
    REPORT10BETABOXPLOT("Report_10", ioi_ch, ord_ioi_ch, metadata_ch, report_ten_ch, unifrac_pairwise_ch)
    REPORT11UPGMA( "Report_11", ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, ioi_ch, ord_ioi_ch, report_eleven_ch)
    REPORT12PERMANOVA("Report_12", ch_norm_qza_table, new_rooted_tree_ch, ch_tax_qza, metadata_ch, ioi_ch, ord_ioi_ch, ch_core_distance, report_twelve_ch)
    REPORT13LEFSE("Report_13", ch_lefse_images, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi_ch)
    REPORT14CITATIONS("Report_14",report_fourteen_ch)

    if(params.mock){
        // Step 1 filter qza table for only Mock Data
        QIIME2_ONLYMOCK(metadata_ch, ch_raw_qza_table, mock_val_ch, ioi_ch
            ).qza.set { ch_only_mock_qza_table }

        // Test 1 Check Sequence Quality
        // Step 1 filter repseqs for only Mock Data
        QIIME2_FILTERSEQS(ch_only_mock_qza_table, rep_seq_ch
            ).qza.set { ch_only_mock_seq }
        // Step 2 QIIME2 quality-control evaluate-seqs
        // in: refrence seqs (made elsewhere); observed seqs from step 1
        QIIME2_EVALUATE_SEQS(refrence_seq_ch, ch_only_mock_seq)

        // Test 2 Check Quality of Samples with known composition
        // Step 2 QIIME2 quality-control evaluate-composition
        // in: expected qza table (made elsewhere); observed qza table from right before
        QIIME2_EVALUATE_COMPOSITION(ch_only_mock_qza_table, ch_tax_qza, reference_table_ch)
    }
}

    