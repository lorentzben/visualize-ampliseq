#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "${projectDir}/ampliseq_results"
params.metadata = "${projectDir}/metadata.tsv"
params.ioi = "treatment"
params.ordioi = "ordered_item_of_interest.csv"
params.outdir = "results"
params.rare = 0
params.controls = ""

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
         profile : ${workflow.profile}
         """
         .stripIndent()

input_ch = Channel.fromPath(params.input, checkIfExists: true)
metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
ioi_ch = Channel.of(params.ioi)
ord_ioi_ch = Channel.fromPath(params.ordioi)
rare_val_ch = Channel.of(params.rare)
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
if (params.controls) {
    controls_ch = Channel.fromPath(params.controls, checkIfExists:false)
}
contam_script_ch = Channel.fromPath("${projectDir}/r_scripts/contam_script.r")
    

workflow {
    ord_ioi = ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch)
    
    if (params.controls) {
        filtered_table = FILTERNEGATIVECONTROL(input_ch, controls_ch, metadata_ch, contam_script_ch)
        qza_table = TSVTOQZA(FILTERNEGATIVECONTROL.out.filtered_table_biom, metadata_ch)
        

        RAREFACTIONPLOT(input_ch, rare_report_ch, qza_table)
        tax_qza = REFORMATANDQZATAX(input_ch)
        (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza)
        COREMETRICPYTHON(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch)
        QZATOTSV(COREMETRICPYTHON.out.vector)
        REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch, FILTERNEGATIVECONTROL.out.filtered_table_tsv)
        graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
        REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
        REPORT03HEATMAP(input_ch, qza_table, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
        REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
        REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
        REPORT06ORDINATION(qza_table, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
        GENERATERAREFACTIONCURVE(metadata_ch, qza_table, input_ch, count_minmax_ch, rare_val_ch)
        REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
        REPORT08RANKEDABUNDANCE(qza_table,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
        REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
        UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
        GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
        REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
        REPORT11UPGMA( qza_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
        REPORT12PERMANOVA(qza_table, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
        LEFSEFORMAT(ioi_ch, qza_table, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
        lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
        REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
        REPORT14CITATIONS(report_fourteen_ch)
    }
    else{
        RAREFACTIONPLOT(input_ch, rare_report_ch)
        tax_qza = REFORMATANDQZATAX(input_ch)
        (graphlan_biom, table_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza)
        COREMETRICPYTHON(metadata_ch, table_qza, input_ch, count_minmax_ch, rare_val_ch)
        QZATOTSV(COREMETRICPYTHON.out.vector)
        REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch)
        graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
        REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch, report_two_ch, report_two_local_ch)
        REPORT03HEATMAP(input_ch, table_qza, tax_qza, metadata_ch, report_three_ch, ioi_ch, ord_ioi)
        REPORT04ALPHATABLE(QZATOTSV.out.vector, ioi_ch, report_four_ch)
        REPORT05ALPHABOXPLOT(QZATOTSV.out.vector, ioi_ch, ord_ioi, metadata_ch, report_five_ch)
        REPORT06ORDINATION(table_qza, input_ch, ioi_ch, ord_ioi, report_six_ch, tax_qza, metadata_ch, COREMETRICPYTHON.out.pcoa, COREMETRICPYTHON.out.vector)
        GENERATERAREFACTIONCURVE(metadata_ch, table_qza, input_ch, count_minmax_ch, rare_val_ch)
        REPORT07RAREFACTION(ioi_ch,ord_ioi,input_ch, report_seven_ch, GENERATERAREFACTIONCURVE.out.rareVector, metadata_ch)
        REPORT08RANKEDABUNDANCE(table_qza,input_ch, ioi_ch, ord_ioi, report_eight_ch, tax_qza, metadata_ch)
        REPORT09UNIFRACHEATMAP(ioi_ch, ord_ioi, metadata_ch, COREMETRICPYTHON.out.distance, report_nine_ch)
        UNCOMPRESSDIVMATS(COREMETRICPYTHON.out.distance, uncompress_script_ch)
        GENERATEUNIFRAC(COREMETRICPYTHON.out.distance, metadata_ch, ioi_ch)
        REPORT10BETABOXPLOT(ioi_ch,ord_ioi,metadata_ch,input_ch, report_ten_ch, GENERATEUNIFRAC.out.pairwise)
        REPORT11UPGMA( table_qza, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, report_eleven_ch)
        REPORT12PERMANOVA(table_qza, input_ch, ioi_ch, ord_ioi, tax_qza, metadata_ch, COREMETRICPYTHON.out.distance, report_twelve_ch)
        LEFSEFORMAT(ioi_ch, table_qza, input_ch, tax_qza, metadata_ch, qiime_to_lefse_ch)
        lefse_dir = LEFSEANALYSIS(LEFSEFORMAT.out.combos,lefse_analysis_ch, plot_clado_file_ch, plot_res_file_ch)
        REPORT13LEFSE(lefse_dir, report_thirteen_ch, report_thirteen_local_ch, ioi_ch, ord_ioi)
        REPORT14CITATIONS(report_fourteen_ch)
    }
}

process FILTERNEGATIVECONTROL{

    publishDir "${params.outdir}/filtered-table", pattern: "*.qza", mode: "copy"
    publishDir "${params.outdir}/filtered-table", pattern: "*.tsv", mode: "copy"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/decontam:1.2' : 'lorentzb/decontam:1.2' }"

    input:
    path 'results'
    path controls
    path "metadata.tsv"
    path control_script
 

    output:
    path("*.biom"), emit: filtered_table_biom
    path("*.tsv"), emit: filtered_table_tsv
    

    script:

    '''
    #!/usr/bin/env bash

    Rscript contam_script.r

    '''

}

process TSVTOQZA{
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:
    path "table.biom"
    path "metadata.tsv"
 

    output:
    path("*.qza"), emit: table_qza
    

    script:

    '''
    #!/usr/bin/env bash
    
    biom add-metadata -i table.biom -o md-table.biom --observation-metadata-fp metadata.tsv

    qiime tools import \
    --input-path md-table.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza
    '''
}

process RAREFACTIONPLOT{
    publishDir "${params.outdir}/rarefaction", pattern: "*.png", mode: "copy"
    publishDir "${params.outdir}/rarefaction", pattern: "*.csv", mode: "copy"
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/tidyverse:4.2.0' : 'lorentzb/tidyverse:4.2.0' }"

    input:
    path 'results'
    path report
    path table
 

    output:
    path("*.png"), emit: rare_images
    path("*.csv"), emit: rare_tabs
    path("*.html"), emit: rare_report

    script:
    def table = table.name != 'NO_FILE' ? "$table" : ''
    '''
    #!/usr/bin/env bash

    cp $table table.qza

    Rscript -e "rmarkdown::render('rarefaction_report.Rmd', output_file='$PWD/rarefaction_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    '''

}

process COREMETRICPYTHON{
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:

    file 'metadata.tsv'
    file 'table.qza'
    path 'results'
    path 'count_table_minmax_reads.py'
    val rare_val
    

    output:

    path("diversity_core/*_pcoa_results.qza")   , emit: pcoa
    path("diversity_core/*_vector.qza")         , emit: vector
    path("diversity_core/*_distance_matrix.qza"), emit: distance
    path("*rarefaction.txt") , emit: depth

    script:

    """
    #!/usr/bin/env python3

    from qiime2.plugins import diversity
    from qiime2 import Metadata
    from qiime2.plugins import feature_table
    from qiime2 import Artifact
    import pandas as pd
    import sys
    import warnings
    import os

    os.mkdir("diversity_core")

    warnings.filterwarnings('ignore')

    metadata = Metadata.load('metadata.tsv')

    unrarefied_table = Artifact.load('table.qza')

    rooted_tree = Artifact.load('results/qiime2/phylogenetic_tree/rooted-tree.qza')

    # if the default value aka use count_table_minmax_reads
    if $rare_val == 0:

        uncompress_table='results/qiime2/abundance_tables/feature-table.tsv'

        # adapted from count_table_minmax_reads.py @author Daniel Straub
        # collected from nf-core/ampliseq
        # read tsv and skip first two rows
        data = pd.read_csv('results/qiime2/abundance_tables/feature-table.tsv', sep="\t", skiprows=[0, 1], header=None)  # count table

        # drop feature ids
        df = data.drop(data.columns[0], axis=1)

        # make sums
        sums = df.sum()

        # we want minimum values
        mindepth = int(sums.min())

        if mindepth > 10000:
            print("Use the sampling depth of " +str(mindepth)+" for rarefaction")
        elif mindepth < 10000 and mindepth > 5000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is quite small for rarefaction")
        elif mindepth < 5000 and mindepth > 1000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is very small for rarefaction")
        elif mindepth < 1000: 
            print("WARNING The sampling depth of "+str(mindepth)+" seems too small for rarefaction")
        else:
            print("ERROR this shouldn't happen")
            exit(1)

        core = diversity.pipelines.core_metrics_phylogenetic(unrarefied_table, rooted_tree, mindepth, metadata)
        file = open("rarefaction.txt", "w")
        file.write(str(mindepth))
        file.close 
    
    # else if user submits the rarefaction depth they want to use based on rarefaction plot
    else: 
        core = diversity.pipelines.core_metrics_phylogenetic(unrarefied_table, rooted_tree, $rare_val, metadata)
        file = open("rarefaction.txt", "w")
        file.write(str($rare_val))
        file.close 

    Artifact.save(core[0], "diversity_core/rarefied_table")
    Artifact.save(core[1], "diversity_core/faith_pd_vector")
    Artifact.save(core[2], "diversity_core/observed_features_vector")
    Artifact.save(core[3], "diversity_core/shannon_vector")
    Artifact.save(core[4], "diversity_core/evenness_vector")
    Artifact.save(core[5], "diversity_core/unweighted_unifrac_distance_matrix")
    Artifact.save(core[6], "diversity_core/weighted_unifrac_distance_matrix")
    Artifact.save(core[7], "diversity_core/jaccard_distance_matrix")
    Artifact.save(core[8], "diversity_core/bray_curtis_distance_matrix")
    Artifact.save(core[9], "diversity_core/unweighted_unifrac_pcoa_results")
    Artifact.save(core[10], "diversity_core/weighted_unifrac_pcoa_results")
    Artifact.save(core[11], "diversity_core/jaccard_pcoa_results")
    Artifact.save(core[12], "diversity_core/bray_curtis_pcoa_results")
    Artifact.save(core[13], "diversity_core/unweighted_unifrac_emperor")
    Artifact.save(core[14], "diversity_core/weighted_unifrac_emperor")
    Artifact.save(core[15], "diversity_core/jaccard_emperor")
    Artifact.save(core[16], "diversity_core/bray_curtis_emperor")  

    """
}

process QZATOTSV{

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:

    path ( diversity )
    
    output:

    path("*_vector.tsv"), emit: vector

    script:

    """
    #!/usr/bin/env python3

    from qiime2.plugins import diversity
    from qiime2 import Metadata
    from qiime2.plugins import feature_table
    from qiime2 import Artifact
    import pandas as pd
    import sys
    import warnings
    import os

    diversity_names = '$diversity'
    diversity_names = diversity_names.split(' ')

    for item in diversity_names:

        diversity_obj = Artifact.load(item)

        Artifact.export_data(diversity_obj,'.')
        artifact_name = item
        filename = str(artifact_name.split('.')[0]+'.tsv')

        os.rename('alpha-diversity.tsv', filename)
    """
}

process ORDERIOI{

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:
    val ioi
    path metadata
    path ord_ioi

    output:

    file 'order_item_of_interest.csv'

    script:

    """
    #!/usr/bin/env python3
    import pandas as pd

    try:
        # if the order ioi exists then write file out and move to following chunks
        read_ord_ioi = pd.read_table("${ord_ioi}",index_col=0,sep=',')
        pd.DataFrame.to_csv(read_ord_ioi, 'order_item_of_interest.csv', index=False)

    except FileNotFoundError:
        # generate the ordered ioi by sorting it and saves it out

        read_metadata = pd.read_table('${metadata}', index_col=0, sep='\t')

        iois = list(pd.Series.unique(read_metadata['${ioi}']))
        ioisdf = pd.DataFrame(iois[0:])
        ioisdf.columns = ['${ioi}']
        ioisdf = ioisdf.sort_values('${ioi}')
        pd.DataFrame.to_csv(ioisdf, 'order_item_of_interest.csv', index=False)
    """
}

process REPORT01BARPLOT{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/html", pattern: "*/*.png", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/microbiome_analyst:1.1' : 'lorentzb/microbiome_analyst:1.1' }"

    //container 'lorentzb/microbiome_analyst:1.1'

    input:

    path 'results'
    file 'metadata.tsv'
    path report
    file 'item_of_interest.csv'
    path table

    output:

    file "01_report_*.html"
    file "barplots/*"
     
    
    script:
    def table = table.name != 'NO_FILE' ? "$table" : ''
    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('01_report_MbA.Rmd', output_file='$PWD/01_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    #Rscript -e "rmarkdown::render('01_report_MbA.Rmd', output_file='$PWD/01_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REFORMATANDQZATAX{

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"
    
    input:
    path 'results'

    output:
    file "taxonomy.qza"

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time

    tax_tab = pd.read_table('results/dada2/ASV_tax_species.tsv', sep='\t')
    tax_tab.columns = tax_tab.columns.str.replace('ASV_ID', 'Feature ID')

    tax_tab[['Domain','Kingdom']] = tax_tab[['Domain','Kingdom']].fillna(value="Unassigned")

    tax_tab['Kingdom'] = 'D_0__' + tax_tab['Kingdom'].astype(str) + ';'
    tax_tab['Phylum'] = 'D_1__' + tax_tab[~tax_tab['Phylum'].isnull()]["Phylum"].astype(str) + ';'
    tax_tab['Class'] = 'D_2__' + tax_tab[~tax_tab['Class'].isnull()]["Class"].astype(str) + ';'
    tax_tab['Order'] = 'D_3__' + tax_tab[~tax_tab['Order'].isnull()]["Order"].astype(str) + ';'
    tax_tab['Family'] = 'D_4__' + tax_tab[~tax_tab['Family'].isnull()]["Family"].astype(str) + ';'
    tax_tab['Genus'] = 'D_5__' + tax_tab[~tax_tab['Genus'].isnull()]["Genus"].astype(str) + ';'
    tax_tab['Species'] = 'D_6__' + tax_tab[~tax_tab['Species'].isnull()]["Species"].astype(str) + ';'

    tax_tab = tax_tab.fillna('')


    cols = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family','Genus', 'Species']
    tax_tab['Taxon'] = tax_tab[cols].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    
    tax_tab = tax_tab.drop(columns=['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family','Genus', 'Species','sequence'])

    colnames  = ['Feature ID', 'Taxon', 'confidence']

    tax_tab = tax_tab.reindex(columns = colnames)

    tax_tab.to_csv("taxonomy.tsv",sep='\t',index=False)
    
    create_tax_qza_command = "qiime tools import \
    --input-path taxonomy.tsv \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --output-path taxonomy.qza"
    result = subprocess.run([create_tax_qza_command], shell=True)

    """
}

process GENERATEBIOMFORGRAPHLAN{
    //publishDir "${params.outdir}/graphlan", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    //container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file metadata 
    val ioi 
    path 'results'
    file "filter_samples.py" 
    file "taxonomy.qza"

    output:
    
    path "biom_tabs/*" 
    file "feature-table.qza"
    

    //TODO add these labels back in
    //label 'process_medium'

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    ioi_set = set(metadata_table[\"${ioi}\"])
    ioi = '${ioi}'

    subprocess.run(['mkdir phylo_trees'], shell=True)
    subprocess.run(['mkdir biom_tabs'], shell=True)

    create_qza_command = "qiime tools import \
    --input-path results/qiime2/abundance_tables/feature-table.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path feature-table.qza"
    result = subprocess.run([create_qza_command], shell=True)

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi
        
        filter_command = "python3 filter_samples.py -m ${metadata} -i ${ioi} -c "+str(item)
        result = subprocess.run([filter_command], shell=True)

        time.sleep(2)

        # adds taxonomic info needed for plotting
        collapse_command = 'qiime taxa collapse \
        --i-table '+str(item)+'-filtered-table.qza \
        --o-collapsed-table collapse-'+str(item)+'-table.qza \
        --p-level 7 \
        --i-taxonomy taxonomy.qza'

        result = subprocess.run([collapse_command], shell=True)

        # exports artifact so that the next step can collect it
        export_command='qiime tools export \
        --input-path collapse-'+str(item)+'-table.qza \
        --output-path collapse-'+str(item)+'-frequency/'
        
        result = subprocess.run([export_command], shell=True)

        # turns feature table into a human-reable format
        biom_command = 'biom convert -i collapse-'+str(item)+\
        '-frequency/feature-table.biom -o otu-'+str(item)+\
        '-table.tsv --to-tsv --header-key taxonomy'

        result = subprocess.run([biom_command], shell=True)

        # formatting the table so that it is in the correct order
        table = pd.read_table(\"otu-"+str(item)+"-table.tsv\", sep='\t', header=1)
        table = table.drop(columns=['taxonomy'])
        table = table.rename(columns={'#OTU ID':'taxonomy'})
        tax = table.pop('taxonomy')
        insertion_site = len(table.columns)
        table.insert(insertion_site, 'taxonomy', tax)
        table.insert(0, 'OTU_ID', np.arange(len(table)))
        table.to_csv('otu-'+str(item)+'-mod-table.tsv', sep='\t', index=False)

        # human readable table into compressed computer-readble format
        biom_format_command='biom convert -i otu-'+str(item)+ \
        '-mod-table.tsv -o '+str(item)+'-otu-table-mod.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy'

        result = subprocess.run([biom_format_command], shell=True)

        result = subprocess.run(['cp '+str(item)+'-otu-table-mod.biom biom_tabs'],shell=True)
    """

}

process RUNGRAPHLAN{
    publishDir "${params.outdir}/graphlan", pattern: "*/*.png", mode: 'copy'
    publishDir "${params.outdir}/html", pattern: "*/*.png", mode: "copy"
    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/py2_test:2.0' : 'lorentzb/py2_test:2.0' }"
    //container "docker://lorentzb/py2_test:2.0"

    input:

    file metadata 
    val ioi
    file "taxonomy.qza"
    file "graph.sh"
    path "biom_tabs/*"
   
    output:
    path "phylo_trees/*" 
    
    
    //label 'process_low'

    script:
    """
    #!/usr/bin/env python2
    import subprocess
    import csv
    import pandas as pd
    import numpy as np 
    import time
    import os

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    ioi_set = set(metadata_table[\"${ioi}\"])
    ioi = '${ioi}'

    os.system('cp biom_tabs/*-otu-table-mod.biom .')

    os.system('mkdir phylo_trees')

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi

        # Outputs the current ioi so that it can be annotated in the graphlan image
        with open('current.txt', 'w') as file:
            file.write(item)

        # bash script call to handle the steps within a conda python 2.7.17 envionment
        generate_image_command = 'bash graph.sh'
        result = os.system(generate_image_command)

    rename_image = 'cp *_image_graph.png phylo_trees/.'
    result = os.system(rename_image)

    rename_pdf_image = 'cp *_image_pdf_graph.png phylo_trees/.'
    result = os.system(rename_pdf_image)


    """

}


process REPORT02GRAPHLANPHYLOGENETICTREE{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/figures", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_02:2.0' : 'lorentzb/r_02:2.0' }"
    //container 'lorentzb/microbiome_analyst:1.1'

    input:

    path "phylo_trees/*" 
    file 'item_of_interest.csv'
    file "02_report.Rmd"
    file "02_report_local.Rmd"

    output:

    file "02_report_*.html"
     
    script:
    if(workflow.profile.contains('local'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('02_report_local.Rmd', output_file='$PWD/02_report_local_$dt.html', output_format='html_document',clean=TRUE,  knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('02_report_local.Rmd', output_file='$PWD/02_report_local_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else if (workflow.profile.contains('slurm'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.html', output_format='html_document',clean=TRUE,  knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else
        error "I'm not sure which to run, you must use local or slurm profiles"
    
}

process REPORT03HEATMAP{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/heatmap", pattern: "*/*.png", mode: "copy"
    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_03:2.0' : 'lorentzb/r_03:2.0' }"

    input: 

    path 'results'
    file 'feature-table.qza'
    file 'taxonomy.qza'
    file 'metadata.tsv'
    file '03_report.Rmd'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'

    output:

    file "03_report_*.html"
    file "03_report_*.pdf"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    cp -L 03_report.Rmd $PWD/03_report_test.Rmd

    Rscript -e "rmarkdown::render('03_report_test.Rmd', output_file='$PWD/03_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('03_report_test.Rmd', output_file='$PWD/03_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REPORT04ALPHATABLE{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_04:2.0' : 'lorentzb/r_04:2.0' }"


    input: 

    path (vectors)
    file 'item_of_interest.csv'
    file '04_report.Rmd'

    output:

    file "04_report_*.html"
    file "04_report_*.pdf"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    cp -L 04_report.Rmd $PWD/04_report_test.Rmd

    Rscript -e "rmarkdown::render('04_report_test.Rmd', output_file='04_report_$dt.html', output_format='html_document', output_dir='$PWD', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('04_report_test.Rmd', output_file='04_report_$dt.pdf', output_format='pdf_document', output_dir='$PWD', clean=TRUE, knit_root_dir='$PWD')"
    '''

}

process REPORT05ALPHABOXPLOT{
    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_05:2.0' : 'lorentzb/r_05:2.0' }"

    input: 

    path (vectors)
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file 'metadata.tsv'
    file '05_report.Rmd'

    output:

    file "05_report_*.html"
    path "alpha_diversity_boxplot/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('05_report.Rmd', output_file='$PWD/05_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    
    '''
}


process REPORT06ORDINATION{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_06:2.0' : 'lorentzb/r_06:2.0' }"

    input: 

    file 'feature-table.qza'
    path 'results'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file '06_report.Rmd'
    file 'taxonomy.qza'
    file 'metadata.tsv'
    path pcoas
    path vectors

    output:

    file "06_report_*.html"
    file "06_report_*.pdf"
    path "beta_diversity_ordination/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process GENERATERAREFACTIONCURVE{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:
    
    file 'metadata.tsv'
    file 'table.qza'
    path 'results'
    path 'count_table_minmax_reads.py'
    val rare_val

    output:
    path("rarefact/observed_features.csv"), emit: rareVector

    script:

    """
    #!/usr/bin/env python3


    from qiime2.plugins.diversity.visualizers import alpha_rarefaction
    from qiime2.plugins import diversity
    from qiime2 import Metadata
    from qiime2.plugins import feature_table
    from qiime2 import Artifact
    import pandas as pd
    import sys
    import warnings
    import os

    table = Artifact.load('table.qza')
    metadata = Metadata.load('metadata.tsv')
    rooted_tree = Artifact.load('results/qiime2/phylogenetic_tree/rooted-tree.qza')

    # if the default value aka use count_table_minmax_reads
    if $rare_val == 0:

        uncompress_table='results/qiime2/abundance_tables/feature-table.tsv'

        # adapted from count_table_minmax_reads.py @author Daniel Straub
        # collected from nf-core/ampliseq
        # read tsv and skip first two rows
        data = pd.read_csv('results/qiime2/abundance_tables/feature-table.tsv', sep="\t", skiprows=[0, 1], header=None)  # count table

        # drop feature ids
        df = data.drop(data.columns[0], axis=1)

        # make sums
        sums = df.sum()

        # we want minimum values
        mindepth = int(sums.min())
        maxdepth = int(sums.max())

        if mindepth > 10000:
            print("Use the sampling depth of " +str(mindepth)+" for rarefaction")
        elif mindepth < 10000 and maxdepth > 5000: 
            print("WARNING The sampling depth of "+str(mindepth)+" is quite small for rarefaction")
        elif mindepth < 5000 and mindepth > 1000:
            print("WARNING The sampling depth of "+str(mindepth)+" is very small for rarefaction")
        elif mindepth < 1000:
            print("WARNING The sampling depth of "+str(mindepth)+" seems too small for rarefaction")
        else:
            print("ERROR this shouldn't happen")
            exit(1)

        # TODO convert this from bash to python
        #check values
        
        if maxdepth > 75000:
            maxdepth = 75000
        
        if maxdepth > 5000:
            maxsteps=250
        else:
            maxsteps=(maxdepth/20)

        rarefact = alpha_rarefaction(table=table, max_depth=maxdepth, phylogeny=rooted_tree, steps=maxsteps)
        file = open("rarefaction.txt", "w")
        file.write(str(mindepth))
        file.close 
    
    # else if user submits the rarefaction depth they want to use based on rarefaction plot
    else: 
        rarefact = alpha_rarefaction(table=table, max_depth=$rare_val, phylogeny=rooted_tree)
        file = open("rarefaction.txt", "w")
        file.write(str($rare_val))
        file.close 

    
    Artifact.export_data(rarefact.visualization,'rarefact')


    """
}

process REPORT07RAREFACTION{
    
    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_07:2.0' : 'lorentzb/r_07:2.0' }"

    input: 

    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    path 'results'
    file '07_report.Rmd'
    path rarefact
    path 'metadata.tsv'

    output:

    file "07_report_*.html"
    file "07_report_*.pdf"
    path "rarefaction_plots/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REPORT08RANKEDABUNDANCE {
    
    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_08:2.0' : 'lorentzb/r_08:2.0' }"

    input: 

    file 'feature-table.qza'
    path 'results'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file '08_report.Rmd'
    file 'taxonomy.qza'
    file 'metadata.tsv'
    
    output:

    file "08_report_*.html"
    file "08_report_*.pdf"
    path "ranked_abundance_curves/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REPORT09UNIFRACHEATMAP{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_09:2.0' : 'lorentzb/r_09:2.0' }"

    //label 'process_medium'

    input:

    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file "metadata.tsv"
    path distances
    file "09_report.Rmd"

    output:

    file "09_report_*.html"
    file "09_report_*.pdf"
    path "unifrac_heatmaps/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process GENERATEUNIFRAC{

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/automate_16_nf:2.0' : 'lorentzb/automate_16_nf:2.0' }"

    input:

    path distances
    file "metadata.tsv"
    file 'item_of_interest.csv'

    output:

    path("*-pairwise.tsv"), emit: pairwise
     
    script:

    '''
    #!/usr/bin/env bash
    IOI=$(cat item_of_interest.csv)

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    qiime diversity beta-group-significance --i-distance-matrix weighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column $IOI \
    --p-pairwise \
    --o-visualization weighted-unifrac.qzv

    qiime tools export \
    --input-path weighted-unifrac.qzv \
    --output-path weighted-unifrac

    cp weighted-unifrac/raw_data.tsv ./weighted-unifrac-pairwise.tsv

    qiime diversity beta-group-significance --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column $IOI \
    --p-pairwise \
    --o-visualization unweighted-unifrac.qzv

    qiime tools export \
    --input-path unweighted-unifrac.qzv \
    --output-path unweighted-unifrac

    cp  unweighted-unifrac/raw_data.tsv ./unweighted-unifrac-pairwise.tsv
    '''
}

process UNCOMPRESSDIVMATS{
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/qiime2r:1.0' : 'lorentzb/qiime2r:2.0' }"

    input:
    path distances
    file "uncompress_diversity.r"

    output:

    path("*_distance.tsv"), emit: distance

    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    R --no-save < uncompress_diversity.r
    '''

}

process REPORT10BETABOXPLOT{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_10:2.0' : 'lorentzb/r_10:2.0' }"

    input:

    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file "metadata.tsv"
    path "results"
    file "10_report.Rmd"
    path distances

    output:

    file "10_report_*.html"
    file "10_report_*.pdf"
    path "beta_div_boxplots/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''

}

process REPORT11UPGMA{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_11:2.0' : 'lorentzb/r_11:2.0' }"

    input:

    file 'feature-table.qza'
    path 'results'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file 'taxonomy.qza'
    file 'metadata.tsv'
    file '11_report.Rmd'

    output:

    file "11_report_*.html"
    file "11_report_*.pdf"
    path "upgma_plots/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}

process REPORT12PERMANOVA{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_12:2.0' : 'lorentzb/r_12:2.0' }"

    input:

    file 'feature-table.qza'
    path 'results'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'
    file 'taxonomy.qza'
    file 'metadata.tsv'
    path distances
    file '12_report.Rmd'

    output:

    file "12_report_*.html"
    path "anosim/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('12_report.Rmd', output_file='$PWD/12_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
    
    
    '''
}

process LEFSEFORMAT {
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/qiime2lefse:1.0' : 'lorentzb/qiime2lefse:1.0' }"

    input:
    val ioi
    file "feature-table.qza" 
    path "results"
    file "taxonomy.qza"
    file "metadata.tsv"
    file "qiime_to_lefse.R" 
      
    output:
    path "combos/*", emit : combos

    script:
    """
    #!/usr/bin/env bash
    mkdir combos
    Rscript qiime_to_lefse.R ${ioi}
    mv lefse_formatted.txt combos/
    """
}


process LEFSEANALYSIS{
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/py2_env:1.0' : 'lorentzb/py2_env:1.0' }"

    input:

    path "combos/*" 
    file "lefse_analysis.sh"
    file plot_clado 
    file plot_res 

    output:
    path "lefse_images/*"

    script:
    """
    #!/usr/bin/env bash
    mkdir lefse_images
    bash lefse_analysis.sh
    """
}

process REPORT13LEFSE{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_13:2.0' : 'lorentzb/r_13:2.0' }"

    input:

    path "lefse_images/*"
    file '13_report.Rmd'
    file '13_report_local.Rmd'
    file 'item_of_interest.csv'
    file 'order_item_of_interest.csv'

    output:

    file "13_report_*.html"
    //file "13_report_*.pdf"
    
         

    script:
    if(workflow.profile.contains('local'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('13_report_local.Rmd', output_file='$PWD/13_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('13_report_local.Rmd', output_file='$PWD/13_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else if (workflow.profile.contains('slurm'))
        '''
        #!/usr/bin/env bash
   
        dt=$(date '+%d-%m-%Y_%H.%M.%S');
        ls -lRh
        echo $PWD
    
        Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
        #Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
        '''
    else
        error "I'm not sure which to run, you must use local or slurm profiles"
    
}

process REPORT14CITATIONS{

    publishDir "${params.outdir}/html", pattern: "*.html", mode: "copy"
    publishDir "${params.outdir}/pdf", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_14:2.0' : 'lorentzb/r_14:2.0' }"

    input:

    file "14_report.Rmd"

    output:

    file "14_report_*.html"
    file "14_report_*.pdf"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''


}


