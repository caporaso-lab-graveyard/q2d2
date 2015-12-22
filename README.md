# q2d2: QIIME 2 is documentation-driven (a public prototype)

This repository is a publicly accessible prototype of some ideas that we're exploring for QIIME 2. For some information on what's coming in QIIME 2, see [@gregcaporaso's blog post, "Toward QIIME 2"](https://qiime.wordpress.com/2015/10/30/toward-qiime-2/).

**This repository is not QIIME 2 - that doesn't exist yet.** The code in this repository is untested and highly experimental. You should not use this for any real analysis, only for exploring the documentation-driven interface and interactive visualizations. If you're interested in using the underlying methods (e.g., UniFrac, Faith PD, ANCOM), those are available in [scikit-bio 0.4.1](http://scikit-bio.org/) and later.

## Installation (conda-based)

We only support installation of Q2D2 using ``conda``. Other methods may work, but we'll keep these instructions up-to-date. This will install a working development environment for Q2D2.

```bash
conda create -n q2d2 -c https://conda.anaconda.org/biocore python=3.5 scikit-bio jupyter pyyaml
source activate q2d2
pip install https://github.com/rossant/ipymd/archive/master.zip
git clone https://github.com/gregcaporaso/q2d2.git
cd q2d2
pip install -e .
```

## Usage: Graphical interface

### Launch the server

Run the following commands to launch the Q2D2/Jupyter server:

```bash
mkdir my-study
cd my-study
q2d2 serve
```

### Load data files

You'll be prompted to load files in the web browser window that opens. The following table describes the input files.

|  File type  | Example location  | Description  |
|---|:-:|---|
| ``sample_metadata``  |  ``q2d2/example-data/keyboard/sample-md.tsv`` | Per-sample metadata (e.g., a QIIME 1 mapping file). These can be validated with [Keemei](http://keemei.qiime.org/).  |
| ``otu-metadata`` |  ``q2d2/example-data/keyboard/q191/otu-md.tsv`` | Per-otu metadata (e.g., a QIIME 1 *taxonomy map*, such as those resulting from running ``assign_taxonomy.py``). |
| ``tree``  |  ``q2d2/example-data/keyboard/q191/rep-set.tre`` | A rooted phylogenetic tree containing all of the OTU IDs in ``unrarefied_biom``. A tree generated by QIIME 1's ``make_phylogeny.py`` (i.e., FastTree2) **will not** work directly, as that will be an unrooted tree. An approach for rooting an unrooted trees is described in scikit-bio's [#1227](https://github.com/biocore/scikit-bio/issues/1227#issue-121285751). |
| ``unrarefied_biom``  |  ``q2d2/example-data/keyboard/q191/otu-table.tsv`` | An unrarefied BIOM table in "classic" (i.e., tab-separated text format). This can be created with a QIIME 1 workflow, and then converted to tab-separated text with ``biom convert --to-tsv``. You should not include an sample or observation/OTU metadata in this file, as that should be provided in ``otu-metadata``. (This will be updated to support any BIOM-formatted file.) |
| ``rarefied_biom``  |  None provided, can be generated using the *Rarefy BIOM Table* workflow. | A rarefied BIOM table (i.e., where all samples have the same number of counts, and that number is greater than 0). The same format restrictions apply as for ``unrarefied_biom``. |

## Q2D2 plugins

A major goal for QIIME 2 will be to simplify the process of integrating and supporting new methods. To do this, we're prototyping a plugin system that can be used to create QIIME 2 compatible workflows. The QIIME 2 codebase itself will primarily be an API, CLI and support for the framework illustrated in Q2D2. Processing methods, visualizations, statistical analyses, and so on will be available as plugins.

The QIIME development group will develop and distribute some of these plugins, but our goal is for it to be straight-forward for anyone to develop, contribute and distribute workflows. **This will be a transition from the QIIME 1 idea that a method "is in QIIME", to the QIIME 2 idea that a method can be used through QIIME.** This will reduce the development burden on our team which will let us focus on core plugins and the underlying framework. It will also empower microbiome software developers to more easily make their software QIIME-compatible. Taken together, this means that our users get access to the latest methods sooner.

## Adding a Q2D2 workflow as a plugin

If you're interested in contributing a workflow to Q2D2 for testing, you will develop it as a Q2D2 plugin. You'll do this in your own repository for that plugin: that workflow won't live in Q2D2 directly, but users will install it as a plugin.

**We're not ready for a lot of plugins yet, as this system is going to be entirely re-written in early 2016.** It would however be very useful to have a few (i.e., around three) for us to test with. The workflows in your plugin should also be short. Use [this](https://github.com/gregcaporaso/q2d2/blob/master/q2d2/plugins/core_diversity/biom-to-adiv.md) and [this](https://github.com/gregcaporaso/q2d2/blob/master/q2d2/plugins/core_diversity/rarefy-biom.md) as your models. This will ensure that as things change in the Q2D2 to QIIME 2 transition, you won't have to re-write a lot of code. Remember, we're in prototyping stage and we're about to transition to alpha development stage. APIs can and will change.  

### Defining a plugin

A plugin has two pieces:
 * a manifest file (``.yml`` file), which describes the plugin including its name and the workflows that it provides
 * one or more workflows (``.md`` files) with YAML front-matter that describes the workflow, including its name, its inputs, and its outputs

Q2D2 includes one built-in plugin, *Core Diversity Analyses*, and you should use the files in [its plugin directory](https://github.com/gregcaporaso/q2d2/tree/master/q2d2/plugins/core_diversity) to model your own plugins. (This plugin will likely transition out of the core repository when we move to QIIME 2, and instead become part of some distribution of directly supported workflows.)

### Installing a plugin

After building a plugin, you can install it in your Q2D2 environment by changing into your plugin's directory (where the ``.yml`` and ``.md`` file(s) are), and running the command:

```bash
q2d2 plugins add .
```

If you want to install your plugin while you're developing it, you should use the command:

```bash
q2d2 plugins add . --link
```

This will allow you to see updates to your plugin's workflows as you edit its ``.md`` files. This is similar to ``pip install -e``.

If you add a plugin while your server is running, you'll need to restart your server to have access to the new workflows.

### Current plugin limitations

For now, plugins can only contain Python 3 code, and can only use libraries that are dependencies of Q2D2 (listed [here](https://github.com/gregcaporaso/q2d2/blob/master/setup.py#L51)). In the future, it will be possible to write these in [any language that has a Jupyter kernel](https://github.com/ipython/ipython/wiki/IPython-kernels-for-other-languages), and they will be able to specify their own dependencies.

Plugins are also currently limited to taking as input and creating as output the file types described under *Load Data Files* above.

## Usage: Command line interface

This will work, but is legacy functionality and will likely be removed. It doesn't offer any benefits over the command line interface

```bash
q2d2 create-study --study-id my-study --sample-metadata-filepath $PWD/example-data/keyboard/sample-md.tsv --otu-metadata-filepath $PWD/example-data/keyboard/q191/otu-md.tsv --tree-filepath $PWD/example-data/keyboard/q191/rep-set.tre --unrarefied-biom-filepath $PWD/example-data/keyboard/q191/otu-table.tsv
ipython notebook my-study/index.md
```
