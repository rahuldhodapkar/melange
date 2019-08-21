# Loading Data from RSEM

A core task pre-requisite of single cell analysis is
the loading of count or transcript abundance data from
a quantitation software. [RSEM (RNA-Seq by Expectation Maximization)](https://www.ncbi.nlm.nih.gov/pubmed/21816040)
is one such software for transcript abundance estimation from aligned
sequencing data. This vignette shows how to load data from multiple files
from `rsem-calculate-expression` into melange, and then into a [Seurat](https://satijalab.org/seurat/)
object for downstream analysis.

## RSEM output files overview
`rsem-calculate-expression` produces two output files by default: `sample_name.genes.results` and
`sample_name.isoforms.results`. Melange loads data at the genes level, and will require the locations
of the `*.genes.results` files.

In R, create two character vectors, one with the names of the cells, and another with the locations
of the corresponding RSEM genes results files.

    # e.g.
    #   cell.names <- c('cell1', 'cell2', 'cell3')
    cell.names

    # e.g.
    #   cell.rsem.file.locations <- c('~/cell1.genes.results', '~/cell2.genes.results', '~/cell3.genes.results')
    cell.rsem.file.locations   

## Loading RSEM data

Load the melange library

    library(melange)

and call `LoadRSEM` to create a `Melange` object containing collated quantification data. By default,
melange will load transcript abundance data as TPM from RSEM input files.

    rsem.melange <- LoadRSEM(cell.names, cell.rsem.file.locations)

This may take some time to run depending on the number of input files provided. If you would like to save the melange
object for analysis later, you can use R's `saveRDS` and `readRDS` utilities:

    # Save Melange object
    saveRDS(rsem.melange, file = "rsem_melange.rds")

    # Load saved Melange object to R workspace
    readRDS(file = "rsem_melange.rds")

## Analysis in Seurat

After data is loaded into melange you can either create a new Seurat object using the data from melange, or
convert a Melange object directly into a Seurat assay, which can be joined with other data from the same cells
in an integrated analysis workflow.

### Create a Seurat Object

As of Seurat v3:

    seurat.object <- CreateSeuratObject(counts = rsem.melange@data, 
                                        project = "MySeuratProject", assay = "RNA")

### Create a Seurat Assay

As of Seurat v3:

    seurat.assay <- SeuratAssayFromMelange(rsem.melange)

To add to an existing Seurat object:

    existing.seurat.object[['AssayName']] <- seurat.assay

