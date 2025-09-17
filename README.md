# mspepsearchr

The `mspepsearchr` package provides an R interface to the MSPepSearch tool from NIST enabling library searches against mass spectral databases in NIST format. The package offers user-friendly access to all library search algorithms available in the NIST MS Search software for EI and MS/MS spectra of small molecules. Results from multiple searches are automatically parsed and returned as a list of data frames.

Supported algorithms and their corresponding functions:
- Identity EI Normal - `IdentitySearchEiNormal()`
- Identity HighRes NoPrecursor - `IdentitySearchHighRes()`
- Identity MS/MS - `IdentitySearchMsMs()`
- Similarity EI Simple - `SimilaritySearchEiSimple()`
- Similarity EI Neutral Loss - `SimilaritySearchEiNeutralLoss()`
- Similarity EI Hybrid - `SimilaritySearchEiHybrid()`
- Similarity MS/MS in EI - `SimilaritySearchMsMsInEi()`
- Similarity MS/MS Hybrid - `SimilaritySearchMsmsHybrid()`

The MSPepSearch tool is included in this package as executables, in compliance with NIST’s distribution policy. While it is uncommon for R packages to contain executables, this approach offers two key advantages: improved convenience for end users and the ability to implement automated testing.

The package is available on Windows, Linux and macOS (on Linux and macOS Wine must be installed and accessible via `$PATH` in order to run the Windows executables).


## Installation

``` r
# Install 'mssearchr' from GitHub:
library(devtools)
install_github("https://github.com/AndreySamokhin/mspepsearchr")
```


## Known Discrepancies Between MS Search and MSPepSearch

Users may notice subtle differencies in results when using MSPepSearch versus MS Search:
- The returned probability (`prob`) may differ slightly (typically less than 0.15%).
- Results from the 'Similarity MS/MS in EI' algorithm (`SimilaritySearchMsMsInEi()`) can differ.
- When using high-resolution mass spectra, results of algorithms originally developed for low-resolution mass spectra (such as 'Identity EI Normal' and 'Similarity EI Simple') may differ.
- RI values are not exported when using the 'Similarity EI Hybrid' algorithm.
- The 'Best matching only' option functions differently: for MSPepSearch, it only works with CAS numbers, whereas MS Search can use compound names if CAS numbers are unavailable.


