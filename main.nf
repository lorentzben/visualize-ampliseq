#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "${projectDir}/ampliseq_results"
params.metadata = "${projectDir}/metadata.tsv"
params.ioi = "treatment"
params.outdir = "results"

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input }
         metadata : ${params.metadata}
         item of interest : ${params.ioi}
         outdir   : ${params.outdir}
         """
         .stripIndent()

input_ch = Channel.fromPath(params.input, checkIfExists: true)
metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
ioi_ch = Channel.of(params.ioi)
report_one_ch = Channel.fromPath("${projectDir}/report_gen_files/01_report_MbA.Rmd")
filter_samples_ch = Channel.fromPath("${projectDir}/python_scripts/filter_samples.py")
graph_sh_ch = Channel.fromPath("${projectDir}/bash_scripts/graph.sh")

workflow {
    REPORT01BARPLOT(input_ch, metadata_ch, report_one_ch, ioi_ch)
    tax_qza = REFORMATANDQZATAX(input_ch)
    (graphlan_biom, taxonomy_qza) = GENERATEBIOMFORGRAPHLAN(metadata_ch, ioi_ch, input_ch, filter_samples_ch, tax_qza)
    graphlan_dir = RUNGRAPHLAN(metadata_ch, ioi_ch, tax_qza, graph_sh_ch, graphlan_biom)
    REPORT02GRAPHLANPHYLOGENETICTREE(graphlan_dir, ioi_ch)
}

process REPORT01BARPLOT{

    publishDir "results/html", pattern: "*.html", mode: "copy"
    publishDir "results/figures", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/microbiome_analyst:1.1' : 'lorentzb/microbiome_analyst:1.1' }"

    //container 'lorentzb/microbiome_analyst:1.1'

    input:

    path 'results'
    path metadata
    path report
    file 'item_of_interest.csv'

    output:

    file "01_report_*.html"
    file "barplots/*"
     
    script:

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
    publishDir "${params.outdir}/graphlan", mode: 'copy'

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

    publishDir "results/html", pattern: "*.html", mode: "copy"
    publishDir "results/figures", pattern: "*/*.png", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/r_02:2.0' : 'lorentzb/r_02:2.0' }"
    //container 'lorentzb/microbiome_analyst:1.1'

    input:

    path "phylo_trees/*" 
    file 'item_of_interest.csv'

    output:

    file "02_report_*.html"
    file "graphlan_trees/*"
     
    script:

    '''
    #!/usr/bin/env bash
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    #Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''
}


process LefseFormat {
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "docker://lorentzb/qiime2lefse:1.0"

    input:
    val ioi from ch_ioi_lefse
    file "table-dada2.qza" from ch_table_lefse_graphlan
    file "rarefied_table.qza" from ch_table_lefse
    file "rooted-tree.qza" from ch_tree_lefse
    file "taxonomy.qza" from ch_tax_lefse
    file metadata from ch_metadata_lefse
    file "qiime_to_lefse.R" from ch_lefse_format_script
      

    output:
    path "combos/*" into ch_paired_lefse_format
    file "table-dada2.qza" into ch_table_report_raw
    file "taxonomy.qza" into ch_tax_report
    file "metadata.tsv" into ( ch_metadata_report, ch_metadata_r01, ch_metadata_r02, ch_metadata_r03, 
    ch_metadata_r04, ch_metadata_r05, ch_metadata_r06, ch_metadata_r07, ch_metadata_r08, ch_metadata_r09, 
    ch_metadata_r10, ch_metadata_r11, ch_metadata_r12, ch_metadata_r13 )

    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash
    mkdir combos
    cp ${metadata} "metadata.tsv"
    Rscript qiime_to_lefse.R ${ioi}
    mv lefse_formatted.txt combos/
    """
}


process LefseAnalysis{
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "docker://lorentzb/py2_env:1.0"

    input:
    path "combos/*" from ch_paired_lefse_format
    file "lefse_analysis.sh" from ch_lefse_analysis_script
    file plot_clado from ch_clado_file
    file plot_res from ch_plot_res

    output:
    path "result/*" into ( ch_lefse_results, ch_lefse_r13 ) 

    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash
    mkdir result
    bash lefse_analysis.sh
    """
}

