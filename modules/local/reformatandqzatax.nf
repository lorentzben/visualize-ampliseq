process REFORMATANDQZATAX{

    tag "$asv_tsv"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    'docker://lorentzb/automate_16_nf:2.0' : 
    'lorentzb/automate_16_nf:2.0' }"
    
    input:
    path (asv_tsv)

    output:
    path("taxonomy.qza"), emit: tax_qza

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time

    tax_tab = pd.read_table(, sep='\t')
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