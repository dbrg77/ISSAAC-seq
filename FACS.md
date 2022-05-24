# ISSAAC-seq by FACS
This section contains information about ISSAAC-seq using a plate-based workflow. FACS is commonly used to sort single nuclei into wells, but handpick should also work.

## About the libraries and files

|          |     Read     |                                         Description                                        |
|:--------:|:------------:|:------------------------------------------------------------------------------------------:|
| __ATAC__ |  Read 1 (R1) |                                          gDNA read                                         |
|          |  Read 2 (R2) |                                          gDNA read                                         |
|          | Index 1 (I1) |                                          i7 index                                          |
|          | Index 2 (I2) |                                          i5 index                                          |
|----------|--------------|--------------------------------------------------------------------------------------------|
|  __RNA__ |  Read 1 (R1) |            The first 10 bp are UMIs, the rest are ignored as they are mostly dT            |
|          |  Read 2 (R2) | cDNA sequence, might have some adaptor contamination depending on how long you sequence it |
|          | Index 1 (I1) |                                          i7 index                                          |
|          | Index 2 (I2) |                                          i5 index                                          |

Check [this page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html#FACS) to see how a step-by-step guide of how the libraries are generated. In this workflow, single nuclei are sorted into invidividual wells containing index primers. The library preparation is done individually. In this case, the well barcode is the cell barcode. The combination of `Index 1 (i7)` and `Index 2 (i5)` defines a cell. There are two common scenarios: 

1. You equence your libraries from a core facility. In this case, you probably need to provide the index information `Index 1 (i7) + Index 2 (i5)`, and the core will demultiplex for you. In this case, you have two fastq files for each cell per modality: `R1.fastq.gz` and `R2.fastq.gz`.

2. You sequence by yourself, and run `bcl2fastq` from `Illumina` to generate fastq files. You can certainly put the `Index 1 (i7) + Index 2 (i5)` information in the `SampleSheet.csv`, and each cell will get demultiplexed by `bcl2fastq`. A simpler way is simply put `NNNNNNNN + NNNNNNNN` in the `SampleSheet.csv`. Then run `bcl2fastq` like this:

```
bcl2fastq -r 4 -p 4 -w 4 --create-fastq-for-index-reads --no-lane-splitting -o /path/to/output/dir
```

Due to a known bug in `bcl2fastq`, the program look for the literal `N` in the index. Therefore, you will get the following four files:

```
Undetermined_S0_L001_R1_001.fastq.gz
Undetermined_S0_L001_R2_001.fastq.gz
Undetermined_S0_L001_I1_001.fastq.gz
Undetermined_S0_L001_I2_001.fastq.gz
```

You can start from there, and use the `Snakefile` to process them.