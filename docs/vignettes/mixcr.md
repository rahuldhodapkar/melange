# Loading Data from MiXCR

[MiXCR](https://mixcr.readthedocs.io/en/master/) is a simple tool for analysis of T- and B- 
cell receptor clones within sequencing
data sets. This data can be useful to identify certain reactive subsets of cells across
different tissue samples, or to identify expanded subpopulations of cells.
This vignette will walk through loading data from MiXCR output files into a format
that can be analyzed alongside other single cell RNA-seq data.

## MiXCR output files overview

MiXCR assembles immune receptor sequences from raw reads and outputs different report
files for each immune receptor chain interrogated. For T-cells, this may correspond to
alpha, beta, gamma, and delta chains. For B- cells, this means the kappa and lambda
light chains.

Melange is designed to load only a single chain type at a time. If you want to
analyze chain pairings, simply load **two** Melange objects, with the same cell
IDs, and the corresponding matrices can be joined.

## Prepare MiXCR files

In R, create two character vectors, one with the names of the cells, and another with the locations
of the corresponding MiXCR output files. Note that you will want to only include files from a
specific clonotype.

    # e.g.
    #   cell.names <- c('cell1', 'cell2', 'cell3')
    cell.names

    # e.g.
    #   cell.rsem.file.locations <- c('~/cell1_clones_TRB.txt', '~/cell2_clones_TRB.txt', '~/cell3_clones_TRB.txt')
    cell.mixcr.file.locations   

## Loading MiXCR files

Load the melange library

    library(melange)

and call `LoadMiXCR` to create a `Melange` object containing the within-cell call frequencies
for each clone. As there may be ambiguity within a single cell as to the specific chain sequence,
multiple candidates may be present in the final matrix per-cell. The rows of the output data
matrix correspond to nucleotide sequences for the CDR3 regions of the receptor output.

    mixcr.melange <- LoadMiXCR(cell.names, cell.rsem.file.locations)

Additional information about each receptor is included in the `meta.features` slot of the
`Melange` object. 

