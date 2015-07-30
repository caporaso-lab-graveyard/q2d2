

## Installation (conda-based)

```bash
conda create -n q2d2 python=3.4 scikit-bio
source activate q2d2
pip install -e ~/code/q2d2
```

Then, follow *Step 3* of the [ipymd install instructions](https://github.com/rossant/ipymd/blob/master/README.md).

## Usage

```bash
q2d2 seqs_to_biom --sequences-filepath $HOME/code/q2d2/example-data/keyboard/forensic-seqs.fna --analysis-root my-analysis
q2d2 biom_to_pcoa --analysis-root my-analysis --metadata-filepath $HOME/code/q2d2/example-data/keyboard/forensic-map.txt --color-by Subject
q2d2 start_server --analysis-root my-analysis
```
