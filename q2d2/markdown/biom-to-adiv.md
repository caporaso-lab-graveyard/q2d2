First, we'll load the table and the sample metadata.

```python
>>> %matplotlib inline
>>> import pandas as pd
>>> from q2d2 import load_sample_metadata
>>> sample_md = load_sample_metadata()
>>> from q2d2 import load_rarefied_table
>>> table = load_rarefied_table()
>>> from q2d2 import load_tree
>>> tree = load_tree()
```

Next, we'll compute alpha diversity for all of the samples in the rarefied BIOM table. By default, we compute Faith PD, Chao1, and observed OTUs. The [scikit-bio alpha diversity API documentation](http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html#module-skbio.diversity.alpha) will allow you to explore the implementation of these methods, and the theory behind the methods are discussed in the [Studying Biological Diversity chapter of *An Introduction to Applied Bioinformatics*](http://readiab.org/book/latest/3/1#3) and in [Lesson 8 of *Gut Check*](https://www.coursera.org/learn/microbiome/). **This step might take a few minutes to run.**

```python
>>> from q2d2 import compute_alphas
>>> alpha_diversities = compute_alphas(table, tree)
```

Next, we can interactively explore the relative alpha diversities of these samples. Choose a  category to plot, optionally add some filtering criteria and a hue to color boxes by, and the metric you'd like to plot. If you'd like to compare multiple metrics, simply make a copy of the following cell. 

```python
>>> from q2d2.wui import interactive_plot_alpha_diversity
>>> interactive_plot_alpha_diversity(sample_md, alpha_diversities)
```
