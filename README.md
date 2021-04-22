---
title: "README.md"
author: "YihanXiao"
date: "4/19/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Pipeline for spatial sequencing data

Pipeline-for-spatial-data is a protocol for processing spatial sequencing data. Currently, it supports [10x-genomics](https://www.10xgenomics.com/products/spatial-gene-expression), allowing both raw data and processed matrixs as input. The flow could also take Processed matrixs from [FISH-seq](https://www.nature.com/articles/s41586-019-1049-y#Sec1), [STARmap](https://science.sciencemag.org/content/361/6400/eaat5691), [Slide-seq](https://science.sciencemag.org/content/363/6434/1463) and [Dbit-seq](https://www.cell.com/cell/pdf/S0092-8674(20)31390-8.pdf) as input files. Processed matrixs consist of count matrix and matrix with location information for each spots obtained from image processing.

The pipeline is composed of for main parts: **Preprocessing**, **Spatial DE genes**, **Spatial Clustering**and **Cell-cell interaction**. `Preprocessing` part includes steps of initialization, visualization, QC, normalization, filter, dimension reduction, etc. `DE genes` part is for the genes with strong difference in spatial expression pattern, in which methods including `SpatialDE`, `SPARK`, `binspect` and `Trendsceek` are collected. `spatial clustering` part investigates sub-domains for tissue heterogeneity through `HMRF`. `Cell-cell interaction` part aim to study the potential gene pair relation and pathways among a population of neighbor cells in space.

## System requirements
* Linux/Unix
* Python (>= 3.0) for MAESTRO snakemake workflow
* R (>= 3.6.1) for MAESTRO R package

## Environment configuration

``` {bash}
cd Environment
Environment.sh
```

## Usage
``` {bash2}
Integrate.sh <countmatrix> <locationfile> <resultfilename> <Demehtod (Trendsceek, Binspect)>
```

##Arguments
Arguments  |  Description
---------  |  -----------
`$1 <countmatrix>` | A count matrix, genes by spots.
`$2 <locationfile>` | A numeric matrix for the spot-positions detected, 2 columns at least. 
`$3 <resultfilename>` | The dir for outputs.
`$4 <DEmethod>` | The method for spatial DE genes.

