

## Installation (conda-based)

```bash
conda create -n q2d2 python=3.4
source activate q2d2
conda install --file ~/code/scikit-bio/ci/conda_requirements.txt
pip install -e ~/code/q2d2
```

## Usage

```bash

# q2d2 seqs-to-biom <analysis-name> <fasta|fastq|fasta.gz|fastq.gz> <exact|split|all|rand|last>
./scripts/q2d2 seqs-to-biom example-analysis example-data/keyboard/forensic-seqs.fna exact
```

Execute the Jupyter notebook that is opened.

```bash

# q2d2 biom-to-diversity <sample-metadata.tsv> <color-by>
./scripts/q2d2 biom-to-pcoa example-analysis example-data/keyboard/forensic-map.txt Subject
```

Execute the Jupyter notebook that is opened.
