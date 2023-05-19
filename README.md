# visualize-ampliseq

A pipeline built using DSL2 and nextflow to visualize the output from a [nf-core/ampliseq](https://github.com/nf-core/ampliseq) analysis. 

### Example Parameters:

```yaml
# paramfile.yaml

input: "/cycle-4/litter-srs | result directory absolute path from nf-core/ampliseq "
ioi: "Treatment | "item of interest" column located in metadata that comparisions will be based on"
metadata: "/cycle-4/litter-srs/metadata.tsv | absolute filepath to metadata file from nf-core/ampliseq analysis"
outdir: "/cycle-4/litter-srs-nc-mock-16500-with-viz | outdir of this analysis where resultfiles will be published"
ordioi: "/cycle-4/litter/ord_ioi.csv | "ordered item of interest" if you prefer a sepecific ordering control vs treatment"
rare: 16500 | rarefaction depth can be determined by examining rarefaction/srs curve
srs: true | does SRS need to be utilized?
negative: "NC | what is the value of negative controls in your ioi col"
mock: "MOCK | what is the value of positive control/mock community in your ioi col"
refSeq: "/reference-community/atcc_20_ref.qza | reference sequences for your mock community in qza format, see below"
refTab: "/reference-community/expected-taxonomy.qza | reference abundances for mock community in qza format, see below"
```

### To run:

One note I suggest is to put the command below and the paramfile in a git repository and log commits as you change paramters for reproducibility 

```bash

$ ml Nextflow/22.04.5

# To run on a slurm server with singularity and using the paramfile above

$ nextflow run lorentzben/visualize-ampliseq \
    -r 2.0 \
    -with-tower \
    -profile slurm,singularity \
    -params-file paramfile.yaml \
    -latest \
    -resume

# To run on a local machine with docker and using the paramfile above

$ nextflow run lorentzben/visualize-ampliseq \
    -r 2.0 \
    -with-tower \
    -profile local,docker \
    -params-file paramfile.yaml \
    -latest \
    -resume

```

### RefSeq

If you have a directory of sequnences you can use the following command:

```bash
cd /reference-community
cat *.fasta > atcc_20_ref.fasta
head /work/sealab/bjl34716/ade/reference-community/atcc_20_ref.fasta

qiime tools import \
  --input-path /work/sealab/bjl34716/ade/reference-community/atcc_20_ref.fasta \
  --output-path /work/sealab/bjl34716/ade/reference-community/atcc_20_ref.qza \
  --type 'FeatureData[Sequence]'
```

### RefTab

```bash

nano reftab.tsv

Taxonomy        MOCK167 MOCK288 MOCK323 MOCK71
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Pseudomonadales;D_4__Moraxellaceae;D_5__Acinetobacter;D_6__baumannii    0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;D_5__Bacillus;D_6__pacificus        0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Pseudomonadales;D_4__Moraxellaceae;D_5__Acinetobacter;__        0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Actinobacteriota;D_2__Actinobacteria;D_3__Bifidobacteriales;D_4__Bifidobacteriaceae;D_5__Bifidobacterium;D_6__adolescentis   0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Clostridiaceae;D_5__Clostridium sensu stricto 1;D_6__beijerinckii 0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Actinobacteriota;D_2__Actinobacteria;D_3__Propionibacteriales;D_4__Propionibacteriaceae;D_5__Cutibacterium;D_6__acnes        0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Deinococcota;D_2__Deinococci;D_3__Deinococcales;D_4__Deinococcaceae;D_5__Deinococcus;D_6__radiodurans R1     0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Enterococcaceae;D_5__Enterococcus;D_6__faecalis    0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Enterobacterales;D_4__Enterobacteriaceae;D_5__Escherichia-Shigella;D_6__coli    0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Campylobacterota;D_2__Campylobacteria;D_3__Campylobacterales;D_4__Helicobacteraceae;D_5__Helicobacter;D_6__pylori    0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Lactobacillaceae;D_5__Lactobacillus;D_6__gasseri   0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Burkholderiales;D_4__Neisseriaceae;D_5__Neisseria;D_6__meningitidis     0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Bacteroidota;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Porphyromonadaceae;D_5__Porphyromonas;D_6__gingivalis  0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Pseudomonadales;D_4__Pseudomonadaceae;D_5__Pseudomonas;__       0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhodobacterales;D_4__Rhodobacteraceae;D_5__Cereibacter;__       0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Actinobacteriota;D_2__Actinobacteria;D_3__Actinomycetales;D_4__Actinomycetaceae;D_5__Actinomyces;D_6__odontolytica   0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Staphylococcales;D_4__Staphylococcaceae;D_5__Staphylococcus;D_6__aureus 0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Staphylococcales;D_4__Staphylococcaceae;D_5__Staphylococcus;D_6__epidermidis    0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Streptococcus;D_6__agalactiae        0.05    0.05    0.05    0.05
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Streptococcaceae;D_5__Streptococcus;D_6__mutans    0.05    0.05    0.05    0.05
```

Taking in the file reftab.tsv you can generate the qza you need as below. You must ensure that your mock community IDs and names are the same or you will get an issue. 

```bash

biom convert \
  -i reftab.tsv \
  -o expected-taxonomy.biom \
  --table-type="OTU table" \
  --to-json

qiime tools import \
 --type FeatureTable[RelativeFrequency] \
 --input-path expected-taxonomy.biom \
 --input-format BIOMV100Format \
 --output-path /work/sealab/bjl34716/ade/reference-community/expected-taxonomy.qza
```