# melange

There are many tools available to generate data for single cell analytics,
and as the field matures, the need to join these data together increases.
`melange` is a toolkit intended to provide intuitive and efficient data
transformation and manipulation utilities for use with the popular
[Seurat](https://github.com/satijalab/seurat) R package.

## Currently Supporting

- [RSEM](https://deweylab.github.io/RSEM/) transcript quantitation by gene

- [Annovar](http://annovar.openbioinformatics.org/) variant annotation hits
    by gene

## Installation Instructions

To install a development release of `melange`, simply clone this repository

    git clone git@github.com:rahuldhodapkar/melange.git

and use `devtools` utilities. For version `3.5.3` of `R`, the invocation is:

    devtools::install("path/to/melange")

At this point, `melange` will be available for use as a normal R package in
your environment by `library(melange)` when required.

To simplify the installation process, you can also navigate to the cloned
melange directory and use:

    make install

which will run the above installation.

## Developing for melange

If you are working with output of a research script and have developed an
in-house solution for data integration with your single-cell sequencing
workflow, please contact us (rahul.dhodapkar at yale.edu). We would love
to include your work as part of melange so others may benefit from your
efforts!
