

## Installation (conda-based)

```bash
conda create -n q2d2 python=3.4
source activate q2d2
conda install --file ~/code/scikit-bio/ci/conda_requirements.txt
pip install -e ~/code/q2d2
```

## Usage

```bash
q2d2 seqs_to_biom --sequences-filepath /Users/caporaso/code/q2d2/example-data/keyboard/forensic-seqs.fna --analysis-root forensic-test
```

Execute the Jupyter notebook that is opened.

```bash
q2d2 biom_to_pcoa --analysis-root example-analysis --metadata-filepath example-data/keyboard/forensic-map.txt --color-by Subject
```

Execute the Jupyter notebook that is opened.
