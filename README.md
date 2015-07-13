

## Installation (conda-based)

```bash
conda create -n q2d2 python=3.4
source activate q2d2
conda install --file ~/code/scikit-bio/ci/conda_requirements.txt
pip install -e ~/code/q2d2
```

## Usage

```bash

q2d2 seqs-to-biom <analysis-name> <fasta|fastq|fasta.gz|fastq.gz> <exact|split|all|random|last>
```

Execute the Jupyter notebook that is opened.

```bash

q2d2 biom-to-diversity <df.gz> <map.tsv> <sampling-depth>
```

Execute the Jupyter notebook that is opened.
