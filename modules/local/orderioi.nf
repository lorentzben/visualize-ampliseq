process ORDERIOI{
    
    tag "$ioi"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'docker://lorentzb/automate_16_nf:2.0' : 
        'lorentzb/automate_16_nf:2.0' }"

    input:
    val (ioi)
    path (metadata)
    path (ord_ioi)

    output:

    file ('order_item_of_interest.csv'), emit: ordered_ioi

    script:

    """
    #!/usr/bin/env python3
    import pandas as pd

    try:
        # if the order ioi exists then write file out and move to following chunks
        read_ord_ioi = pd.read_table("${ord_ioi}",sep=',')
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