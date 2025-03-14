# STRMPS
The `STRMPS` package is designed to extract and collect the short tandem repeat (STR) information from the `fastq` files produced by massively parallel sequencing (MPS). 

## Installation

The `STRMPS`-package depends on `R` (>= 4.4), `methods`, `utils`, `tidyr`, `tibble`, `dplyr`, `stringr`, `purrr`, `parallel`, as well as the bioconductor packages `Biostrings` (>= 2.74.1), `pwalign` (>= 1.2.0), `ShortRead` (>= 1.64.0), and `IRanges` (>= 2.40.1). 

Version 0.5.8 is available on CRAN, but to get the newest version devtools is needed to install the package from github. From R, run the following commands:  

```r
install.packages("devtools")
devtools::install_github("svilsen/STRMPS")
```

NB: The newest version of `STRMPS` (>= 0.6.8) relies `Biostrings` which has been split into two parts `Biostrings` and `pwalign`, i.e., both packages need to available and up-to-date.

## License

This project is licensed under the MIT License.
