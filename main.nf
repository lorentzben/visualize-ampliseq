#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
params.input = "${projectDir}/ampliseq_results"
params.metadata = "${projectDir}/metadata.tsv"
params.ioi = "treatment"
params.ordioi = "ordered_item_of_interest.csv"
params.outdir = "results"
params.rare = 0
params.controls = ""
params.srs = false
params.mock = ""
params.negative = ""

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input }
         metadata : ${params.metadata}
         item of interest : ${params.ioi}
         ordered item of interest : ${params.ordioi}
         outdir   : ${params.outdir}
         rarefaction depth : ${params.rare}
         controls: ${params.controls}
         srs: ${params.srs}
         mock samples: ${params.mock}
         negative control: ${params.negative}
         profile : ${workflow.profile}
         """
         .stripIndent()

input_ch = Channel.fromPath(params.input, checkIfExists: true)
metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
ioi_ch = Channel.of(params.ioi)
ord_ioi_ch = Channel.fromPath(params.ordioi)
rare_val_ch = Channel.of(params.rare)
mock_val_ch = Channel.of(params.mock)
nc_val_ch = Channel.of(params.negative)
rare_report_ch = Channel.fromPath("${projectDir}/r_scripts/rarefaction_report.Rmd")
report_one_ch = Channel.fromPath("${projectDir}/report_gen_files/01_report_MbA.Rmd")
filter_samples_ch = Channel.fromPath("${projectDir}/python_scripts/filter_samples.py")
graph_sh_ch = Channel.fromPath("${projectDir}/bash_scripts/graph.sh")
report_two_ch = Channel.fromPath("${projectDir}/report_gen_files/02_report.Rmd")
report_two_local_ch = Channel.fromPath("${projectDir}/report_gen_files/02_report_local.Rmd")
report_three_ch = Channel.fromPath("${projectDir}/report_gen_files/03_report.Rmd")
report_four_ch = Channel.fromPath("${projectDir}/report_gen_files/04_report.Rmd")
report_five_ch = Channel.fromPath("${projectDir}/report_gen_files/05_report.Rmd")
count_minmax_ch = Channel.fromPath("${projectDir}/python_scripts/count_table_minmax_reads.py")
report_six_ch = Channel.fromPath("${projectDir}/report_gen_files/06_report.Rmd")
report_six_b_ch = Channel.fromPath("${projectDir}/report_gen_files/06b_report.Rmd")
report_seven_ch = Channel.fromPath("${projectDir}/report_gen_files/07_report.Rmd")
report_eight_ch = Channel.fromPath("${projectDir}/report_gen_files/08_report.Rmd")
report_nine_ch = Channel.fromPath("${projectDir}/report_gen_files/09_report.Rmd")
report_ten_ch = Channel.fromPath("${projectDir}/report_gen_files/10_report.Rmd")
report_eleven_ch = Channel.fromPath("${projectDir}/report_gen_files/11_report.Rmd")
report_twelve_ch = Channel.fromPath("${projectDir}/report_gen_files/12_report.Rmd")
qiime_to_lefse_ch = Channel.fromPath("${projectDir}/r_scripts/qiime_to_lefse.R")
lefse_analysis_ch = Channel.fromPath("${projectDir}/bin/lefse_analysis.sh")
plot_clado_file_ch = Channel.fromPath("${projectDir}/python_scripts/plot_cladogram.py")
plot_res_file_ch = Channel.fromPath("${projectDir}/python_scripts/plot_res.py")
report_thirteen_ch = Channel.fromPath("${projectDir}/report_gen_files/13_report.Rmd")
report_thirteen_local_ch = Channel.fromPath("${projectDir}/report_gen_files/13_report_local.Rmd")
report_fourteen_ch = Channel.fromPath("${projectDir}/report_gen_files/14_report.Rmd")
uncompress_script_ch = Channel.fromPath("${projectDir}/r_scripts/uncompress_diversity.r")
srs_curve_ch = Channel.fromPath("${projectDir}/r_scripts/srs_curve.rmd")
srs_min_max_ch = Channel.fromPath("${projectDir}/python_scripts/my_count_table_min_max.py")
if (params.controls) {
    controls_ch = Channel.fromPath(params.controls, checkIfExists:false)
}
contam_script_ch = Channel.fromPath("${projectDir}/r_scripts/contam_script.r")
raw_tsv_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.tsv")
raw_biom_table_ch = Channel.fromPath(params.input+"/qiime2/abundance_tables/feature-table.biom")
rooted_tree_ch = Channel.fromPath(params.input+"results/qiime2/phylogenetic_tree/rooted-tree.qza")
    
include { TSVTOQZA; TSVTOQZA as TSVTOQZA2 } from "${projectDir}/modules/local/tsvtoqza.nf"
include { QIIME2_FILTERSAMPLES as QIIME2_FILTERNC; QIIME2_FILTERSAMPLES as QIIME2_FILTERMOCK } from "${projectDir}/modules/local/qiime2_filtersamples.nf"
include { QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_NC; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_MOCK; QIIME2_EXPORT_ABSOLUTE as QIIME2_EXPORT_ABSOLUTE_CORE  } from "${projectDir}/modules/local/qiime2_export_absolute.nf"
*/

WorkflowMain.initialise(workflow, params, log)

include { VISUALIZEAMPLISEQ } from "${projectDir}/workflows/visualize-ampliseq.nf"

workflow VISUALIZE_AMPLISEQ{
    VISUALIZEAMPLISEQ()
}

workflow{
    VISUALIZE_AMPLISEQ()




    /*
    TODO import these steps into a workflow
    ord_ioi = ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch)
                
            
        //TODO update these values and reference them for downstream processes
        CLEANUPRAWTSV(raw_tsv_table_ch)
        raw_tsv_table = CLEANUPRAWTSV.out.raw_table_tsv
        raw_mba_table = CLEANUPRAWTSV.out.raw_MbA_table_tsv

        CLEANUPRAWQZA(raw_biom_table_ch)
        raw_qza_table = CLEANUPRAWQZA.out.raw_table_qza

        filtered_qza_table = CLEANUPRAWQZA.out.raw_table_qza
        filtered_tsv_table = CLEANUPRAWTSV.out.raw_table_tsv
            
        if(params.controls){
            filtered_table = FILTERNEGATIVECONTROL(input_ch, raw_tsv_table, controls_ch, metadata_ch, contam_script_ch, nc_val_ch)
            tsv_map_1 = FILTERNEGATIVECONTROL.out.filtered_table_biom.map{
                it ->  [ [id: "Filtered-NC-Biom"], it ]
            }

            TSVTOQZA(tsv_map_1, metadata_ch)

            qza_filt_table = TSVTOQZA.out.qza.map{it.last()}

            QIIME2_FILTERNC(metadata_ch, qza_filt_table, nc_val_ch, ioi_ch)

            filtered_qza_table = QIIME2_FILTERNC.out.qza
            QIIME2_EXPORT_ABSOLUTE_NC(QIIME2_FILTERNC.out.qza)
            filtered_tsv_table = QIIME2_EXPORT_ABSOLUTE_NC.out.tsv
        }

        if(params.mock){
            QIIME2_FILTERMOCK(metadata_ch, filtered_qza_table, mock_val_ch, ioi_ch)
            QIIME2_EXPORT_ABSOLUTE_MOCK(QIIME2_FILTERMOCK.out.qza)
            filtered_qza_table = QIIME2_FILTERMOCK.out.qza
            filtered_tsv_table = QIIME2_EXPORT_ABSOLUTE_MOCK.out.tsv
        }

        if(params.srs){

            //if not NC and not Mock filtered_tsv_table should be same as raw table
            //TODO test this
            SRSCURVE(filtered_qza_table, filtered_tsv_table, srs_curve_ch, srs_min_max_ch)
            //if not NC and not Mock filtered_tsv_table should be same as raw table
            //TODO test this
            SRSNORMALIZE(filtered_tsv_table, SRSCURVE.out.min_val, params.rare)

            tsv_map_2 = SRSNORMALIZE.out.biom_normalized.map{
                it ->  [ [id: "SRS-Normalized-Biom"], it ]
            }

            TSVTOQZA2(tsv_map_2, metadata_ch)

            // what we will want to use for most analysis
            norm_qza_table = TSVTOQZA2.out.qza.map{it.last()}
            norm_tsv_table = SRSNORMALIZE.out.tsv_normalized

        }

        
        tax_qza = REFORMATANDQZATAX(input_ch)

        graphlan_biom = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza, filtered_qza_table, nc_val_ch, mock_val_ch)

            
        COREMETRICPYTHON(metadata_ch, filtered_qza_table, norm_qza_table, rooted_tree_ch, rare_val_ch)
        norm_qza_table = COREMETRICPYTHON.out.rare_table
        QIIME2_EXPORT_ABSOLUTE_CORE(COREMETRICPYTHON.out.rare_table)
        norm_tsv_table = QIIME2_EXPORT_ABSOLUTE_CORE.out.tsv
        */
            
        //TODO update these downstream processes

        /*
        REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, SRSNORMALIZE.out.tsv_normalized, tsv_table)
        graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
        REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
        REPORT03HEATMAP(input_ch, COREMETRICPYTHON.out.rare_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
        REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
        REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
        REPORT06ORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
        REPORT06BNMDSORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_b_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
        //TODO update this curve to SRS Curve
        if (!params.srs){
            GENERATERAREFACTIONCURVE(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch, FILTERNEGATIVECONTROL.out.filtered_table_tsv)
            REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
        }
        REPORT08RANKEDABUNDANCE(COREMETRICPYTHON.out.rare_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
        REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
        UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
        GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
        REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
        REPORT11UPGMA( COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
        REPORT12PERMANOVA(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
        LEFSEFORMAT(ioi_ch, qza_table, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
        lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
        REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
        REPORT14CITATIONS(report_fourteen_ch)
        */

    /*
    if (params.srs){
        if (params.controls) {
            //TODO update these values and reference them for downstream processes
            
            

            filtered_qza_table = []
            filtered_tsv_table = []

            CLEANUPRAWTSV(raw_tsv_table_ch)
            raw_tsv_table = CLEANUPRAWTSV.out.raw_table_tsv
            raw_mba_table = CLEANUPRAWTSV.out.raw_MbA_table_tsv

            CLEANUPRAWQZA(raw_biom_table_ch)
            raw_qza_table = CLEANUPRAWQZA.out.raw_table_qza

            filtered_table = FILTERNEGATIVECONTROL(input_ch, controls_ch, metadata_ch, contam_script_ch)

            tsv_map_1 = FILTERNEGATIVECONTROL.out.filtered_table_biom.map{
                it ->  [ [id: "Filtered-NC-Biom"], it ]
            }

            TSVTOQZA(tsv_map_1, metadata_ch)

            qza_filt_table = TSVTOQZA.out.qza.map{it.last()}

            QIIME2_FILTERNC(metadata_ch, qza_filt_table , nc_val_ch, ioi_ch)

            if(params.mock){
                QIIME2_FILTERMOCK(metadata_ch, QIIME2_FILTERNC.out.qza, mock_val_ch, ioi_ch)
                QIIME2_EXPORT_ABSOLUTE(QIIME2_FILTERMOCK.out.qza)
                qza_filt_table = QIIME2_FILTERMOCK.out.qza
                tsv_table = QIIME2_EXPORT_ABSOLUTE.out.tsv
            } else if(params.negative) {
                qza_filt_table = QIIME2_FILTERNC.out.qza
                QIIME2_EXPORT_ABSOLUTE(QIIME2_FILTERNC.out.qza)
                tsv_table = QIIME2_EXPORT_ABSOLUTE.out.tsv
            }

              
            //TODO filter Negative Control Samples and Mock Community Samples
            //TODO export nc/mock filtered table to TSV
        
            SRSCURVE(qza_filt_table, tsv_table, input_ch, srs_curve_ch, srs_min_max_ch)
            tax_qza = REFORMATANDQZATAX(input_ch)
            (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza, qza_filt_table)
    
            SRSNORMALIZE(tsv_table, input_ch, SRSCURVE.out.min_val, params.rare)
            
            tsv_map_2 = SRSNORMALIZE.out.biom_normalized.map{
                it ->  [ [id: "SRS-Normalized-Biom"], it ]
            }

            TSVTOQZA2(tsv_map_2, metadata_ch)

            norm_qza_table = TSVTOQZA2.out.qza.map{it.last()}
            
            COREMETRICPYTHON(metadata_ch, norm_qza_table, input_ch, count_minmax_ch, rare_val_ch)
            
            //TODO make this a module too?
            QZATOTSV(COREMETRICPYTHON.out.vector)
            
            REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, SRSNORMALIZE.out.tsv_normalized, tsv_table)
            graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
            REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
            REPORT03HEATMAP(input_ch, COREMETRICPYTHON.out.rare_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
            REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
            REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
            REPORT06ORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            REPORT06BNMDSORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_b_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            //TODO update this curve to SRS Curve
            //GENERATERAREFACTIONCURVE(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch, FILTERNEGATIVECONTROL.out.filtered_table_tsv)
            //REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
            REPORT08RANKEDABUNDANCE(COREMETRICPYTHON.out.rare_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
            REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
            UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
            GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
            REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
            REPORT11UPGMA( COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
            REPORT12PERMANOVA(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
            LEFSEFORMAT(ioi_ch, qza_filt_table, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
            lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
            REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
            REPORT14CITATIONS(report_fourteen_ch)
        } else {

            tax_qza = REFORMATANDQZATAX(input_ch)

            (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza, [])

            SRSCURVE(table_qza, [] , input_ch, srs_curve_ch, srs_min_max_ch)
            
            //TODO Filter Mock Community Samples
            //TODO export QZA from last process for tsv and update downstream processes

            //if we pass in empty braces, then the process will read in: "results/qiime2/abundance_tables/feature-table.biom"
            SRSNORMALIZE( [] , input_ch, SRSCURVE.out.min_val, params.rare)
            
            tsv_map = SRSNORMALIZE.out.biom_normalized.map{
                it ->  [ [id: "SRS-Normalized-Biom"], it ]
            }

            TSVTOQZA(tsv_map, metadata_ch)

            norm_qza_table = TSVTOQZA.out.qza.map{it.last()}
           
            COREMETRICPYTHON(metadata_ch, norm_qza_table, input_ch, count_minmax_ch, rare_val_ch)
            
            QZATOTSV(COREMETRICPYTHON.out.vector) 
            //empty [] in place of filtered table from decontam
            REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, [])
            graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
            REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
            REPORT03HEATMAP(input_ch, COREMETRICPYTHON.out.rare_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
            REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
            REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
            REPORT06ORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            REPORT06BNMDSORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_b_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            
            REPORT08RANKEDABUNDANCE(COREMETRICPYTHON.out.rare_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
            REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
            UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
            GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
            REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
            REPORT11UPGMA( COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
            REPORT12PERMANOVA(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
            LEFSEFORMAT(ioi_ch, table_qza, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
            lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
            REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
            REPORT14CITATIONS(report_fourteen_ch)
        }

    } else {
        if (params.controls) {
            filtered_table = FILTERNEGATIVECONTROL(input_ch, controls_ch, metadata_ch, contam_script_ch)

            tsv_map_1 = FILTERNEGATIVECONTROL.out.filtered_table_biom.map{
                it ->  [ [id: "Filtered-NC-Biom"], it ]
            }

            TSVTOQZA(tsv_map_1, metadata_ch)
            
            qza_table = TSVTOQZA.out.qza.map{it.last()}

            //TODO if params.mock filter out those and control samples turn to TSV and update downstream
           
            RAREFACTIONPLOT(input_ch, rare_report_ch, qza_table)
            tax_qza = REFORMATANDQZATAX(input_ch)
            (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza, qza_table)
            COREMETRICPYTHON(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch)
            QZATOTSV(COREMETRICPYTHON.out.vector)
            REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, FILTERNEGATIVECONTROL.out.filtered_table_tsv)
            graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
            REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
            REPORT03HEATMAP(input_ch, COREMETRICPYTHON.out.rare_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
            REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
            REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
            REPORT06ORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            REPORT06BNMDSORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_b_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            GENERATERAREFACTIONCURVE(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch, FILTERNEGATIVECONTROL.out.filtered_table_tsv)
            REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
            REPORT08RANKEDABUNDANCE(COREMETRICPYTHON.out.rare_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
            REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
            UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
            GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
            REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
            REPORT11UPGMA( COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
            REPORT12PERMANOVA(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
            LEFSEFORMAT(ioi_ch, qza_table, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
            lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
            REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
            REPORT14CITATIONS(report_fourteen_ch)
        }
        else{
            //TODO if mock remove and export to TSV and update downstream
            empty_table = ord_ioi_ch
            RAREFACTIONPLOT(input_ch, rare_report_ch, [])
            tax_qza = REFORMATANDQZATAX(input_ch)
            (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza, [])
            COREMETRICPYTHON(metadata_ch, table_qza, input_ch, count_minmax_ch, rare_val_ch)
            QZATOTSV(COREMETRICPYTHON.out.vector)
            REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, [])
            graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
            REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
            REPORT03HEATMAP(input_ch, COREMETRICPYTHON.out.rare_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
            REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
            REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
            REPORT06ORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            REPORT06BNMDSORDINATION(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, report_six_b_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
            GENERATERAREFACTIONCURVE(metadata_ch, table_qza, input_ch, count_minmax_ch, rare_val_ch, [])
            REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
            REPORT08RANKEDABUNDANCE(COREMETRICPYTHON.out.rare_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
            REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
            UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
            GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
            REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
            REPORT11UPGMA( COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
            REPORT12PERMANOVA(COREMETRICPYTHON.out.rare_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
            LEFSEFORMAT(ioi_ch, table_qza, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
            lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
            REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
            REPORT14CITATIONS(report_fourteen_ch)
        }

    }
    */
}
