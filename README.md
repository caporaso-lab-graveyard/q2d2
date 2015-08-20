# q2d2: QIIME 2 is documentation-driven (a public prototype)

This repository is a publicly accessible prototype of some ideas that we're exploring for QIIME 2. **This repository is not QIIME 2 - that doesn't exist yet.** The code in this repository is untested and highly experimental. You should not use this for any real analysis, only for exploring the documentation-driven interface.

## Installation (conda-based)

```bash
conda create -n q2d2 python=3.4 scikit-bio ipython-notebook
source activate q2d2
git clone https://github.com/gregcaporaso/q2d2.git
cd q2d2
pip install -e .
```

## Usage

Run the following commands to generate analysis notebooks:

```bash
q2d2 seqs_to_biom --sequences-filepath $PWD/example-data/keyboard/seqs.fna --analysis-root my-analysis
q2d2 rarefy_biom --analysis-root my-analysis
q2d2 biom_to_pcoa --analysis-root my-analysis --metadata-filepath $PWD/example-data/keyboard/sample-md.tsv --color-by Subject
q2d2 biom_to_adiv --analysis-root my-analysis --metadata-filepath $PWD/example-data/keyboard/sample-md.tsv --collated-alpha-filepath $PWD/example-data/keyboard/q191/faith-pd-collated.tsv
q2d2 biom_to_taxa_plots --analysis-root my-analysis --metadata-filepath $PWD/example-data/keyboard/sample-md.tsv --otu-metadata-filepath $PWD/example-data/keyboard/q191/otu-md.tsv --otu-table-filepath $PWD/example-data/keyboard/q191/otu-table.tsv
q2d2 start_server --analysis-root my-analysis
ipython notebook my-analysis/index.md
```
