# Seurat v5 Cross-Species Ortholog Converter

## Overview
This R script provides a robust, highly optimized pipeline to convert single-cell RNA sequencing (scRNA-seq) count matrices between species (e.g., Zebrafish to Human) directly within **Seurat v5** objects. 

It was specifically engineered to overcome two major bottlenecks in cross-species single-cell analysis:
1. **Ensembl Server Firewalls:** It replaces `biomaRt` with `gprofiler2` to completely bypass Ensembl Web Application Firewall (WAF) 403 Forbidden errors and 60-second timeouts.
2. **Seurat v5 `LogMap` Strictness:** It utilizes high-speed sparse matrix multiplication to mathematically aggregate duplicate orthologs (paralogs), preventing the fatal `invalid class "LogMap" object: Duplicate rownames not allowed` error.

## Key Features
* **High-Performance Aggregation:** Replaces buggy, slow `for`-loop row deletions with vectorized `Matrix::fac2sparse` cross-multiplication, instantly summing duplicate gene counts across tens of thousands of cells without index-shifting bugs.
* **API Stability:** Uses the g:Profiler server network (`gorth`) for stable, firewall-free ortholog mapping.
* **Non-Destructive Integration:** Appends the translated count matrix as a distinct, new assay (`RNA_converted`) rather than overwriting the original. This preserves the object's structural integrity, alongside all prior metadata, PCA, and UMAP reductions.
* **Batch Processing:** Supports iterating over a vector of multiple `.rds` files automatically.

## Prerequisites
Ensure the following R packages are installed in your environment. The script is optimized for Seurat v5.

```R
install.packages(c("dplyr", "ggplot2", "xlsx", "gprofiler2", "Matrix", "progress", "future"))
setRepositories(ind = 1:3) # Make sure BiocManager is available if needed for dependencies
install.packages("Seurat")
```

## Usage
1. Open the script and locate the `USER INPUT` section at the top.
2. Define the list of Seurat `.rds` files you wish to convert. 
3. Set your target species using the standard g:Profiler organism code (first letter of genus + species name, all lowercase).

```R
##########################################USER INPUT########################################
# Provide the list of RDS file names here.
file.list <- c("Neurons_zebrafish_ES_complete_081926.rds")

# Target species code (e.g., hsapiens for Human, mmusculus for Mouse)
species.name.to.convert <- "hsapiens" 
############################################################################################
```

4. Source or run the script. The console will output progress bars tracking the matrix conversion and the sparse matrix aggregation steps.

## Expected Output
The script automatically generates and saves a new RDS file in your working directory. The filename will be prefixed with `Converted_`, the target species code, and the original filename.

*Example Output File:* `Converted_hsapiens_Neurons_zebrafish_ES_complete_081926.rds`

Inside the new Seurat object, you will find a new assay named `RNA_converted` containing the translated counts. You can set this as your active assay for downstream integration, comparative analysis, or cross-species label transfer:

```R
DefaultAssay(zebra) <- "RNA_converted"
```

## Troubleshooting: The Dataframe `NULL` Error

If the script crashes during the `data.frame()` construction step with an error resembling:
`arguments imply differing number of rows: 1, 0, 19054`

**Cause:** 
Different versions of the `gprofiler2` package occasionally change the exact column name used to store the input genes. The script defaults to extracting `ortho_raw$input`. If your version of the package returns `NULL` for that column, the dataframe fails to build because it is missing the input vector.

**Solution:**
Open the script, locate the `data.frame()` block inside the `count_matrix_converter` function, and change `$input` to `$initial_alias`.

```R
# Change this:
zebrafish_gene = ortho_raw$input,

# To this:
zebrafish_gene = ortho_raw$initial_alias,
```

## Common Species Codes
| Target Species | Input Code |
| :--- | :--- |
| **Human** | `hsapiens` |
| **Mouse** | `mmusculus` |
| **Zebrafish** | `drerio` |
| **Rat** | `rnorvegicus` |
| **Macaque** | `mmulatta` |

To view the complete, up-to-date registry of all supported species, run `names(gprofiler2::get_version_info()$sources)` in your R console.
