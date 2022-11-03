# visualize-ampliseq

A pipeline built using DSL2 and nextflow to visualize the output from a [nf-core/ampliseq](https://github.com/nf-core/ampliseq) analysis. 

To run:
```bash
$ nextflow run lorentzben/visualize-ampliseq -r dev --input "cec_raw_nf" --ioi "condition" --metadata "cec_raw_nf/all_days_sbm_cec_nf_treatment_metadata.tsv" -profile local,docker

$ nextflow run lorentzben/visualize-ampliseq -r dev --input "cec_raw_nf" --ioi "condition" --metadata "cec_raw_nf/all_days_sbm_cec_q2_metadata.tsv" -profile slurm,singularity

$ ml Nextflow/22.04.5
```