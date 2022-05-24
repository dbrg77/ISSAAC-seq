# ISSAAC-seq by Droplet
This section contains information about ISSAAC-seq using the droplet workflow. Here, the [10x Genomics Single Cell ATAC](https://www.10xgenomics.com/products/single-cell-atac) kit was used, but any droplet systems with a Nextera capture sequence should work, such as [Bio-Rad](https://www.bio-rad.com/en-us/product/surecell-atac-seq-library-prep-kit?ID=PEXSR1MC1ORV) and [HyDrop](https://hydrop.aertslab.org/).

## About the library and files

|          |     Read     |                                         Description                                        |
|:--------:|:------------:|:------------------------------------------------------------------------------------------:|
| __ATAC__ |  Read 1 (R1) |                                          gDNA read                                         |
|          |  Read 2 (R2) |                                          gDNA read                                         |
|          | Index 1 (I1) |                           Sample index (not needed for analysis)                           |
|          | Index 2 (I2) |                                  i5 index (cell barcodes)                                  |
|----------|--------------|--------------------------------------------------------------------------------------------|
|  __RNA__ |  Read 1 (R1) | cDNA sequence, might have some adaptor contamination depending on how long you sequence it |
|          |  Read 2 (R2) |            The first 10 bp are UMIs, the rest are ignored as they are mostly dT            |
|          | Index 1 (I1) |                           Sample index (not needed for analysis)                           |
|          | Index 2 (I2) |                                  i5 index (cell barcodes)                                  |

Check [this page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html#Droplet) to see how a step-by-step guide of how the libraries are generated. In this workflow, single nuclei are caputred using the droplet microfluidics after in situ treatment. Since we use the `10x Genomics Single Cell ATAC` kit for both RNA and ATAC libraries, the sequencing configuration is the same for both modalities:

```
> 50 cycles for Read 1 (R1)
> 50 cycles for Read 2 (R2)
8-10 cycles for I1 (i7) <-- This is the sample index
16 cycles for I2 (i5) <-- This is the cell barcode
```

That means you sequence the libraries as if they are 10x scATAC-seq libraries. Only `R1`, `R2` and `I2` are needed for the analysis. You can basically process them using the `Snakefile` provided in this repository.