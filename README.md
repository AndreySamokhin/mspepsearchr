# mspepsearchr

The `mspepsearchr` package provides an R interface to the MSPepSearch tool from
NIST enabling library searches against mass spectral databases in NIST format.
The package offers user-friendly access to all library search algorithms
available in the NIST MS Search software for EI and MS/MS spectra of small
molecules. Results from multiple searches are automatically parsed and returned
as a list of data frames.

Supported algorithms and corresponding functions:
* Identity EI Normal - `IdentitySearchEiNormal()`
* Identity HighRes NoPrecursor - `IdentitySearchHighRes()`
* Identity MS/MS - `IdentitySearchMsMs()`
* Similarity EI Simple - `SimilaritySearchEiSimple()`
* Similarity EI Neutral Loss - `SimilaritySearchEiNeutralLoss()`
* Similarity EI Hybrid - `SimilaritySearchEiHybrid()`
* Similarity MS/MS in EI - `SimilaritySearchMsMsInEi()`
* Similarity MS/MS Hybrid - `SimilaritySearchMsmsHybrid()`

The MSPepSearch tool is included in this package as executable files, in
compliance with NIST's distribution policy. Although it is uncommon for R
packages to include executables, this approach provides two key advantages:
improved convenience for end users and the ability to perform automated testing.

The package supports external process parallelization via the `parallel`
package, where each R worker runs an independent instance of MSPepSearch.

The package is available on Windows, Linux and macOS (on Linux and macOS Wine
must be installed and accessible via `$PATH` in order to run the Windows
executables).


## Installation

``` r
# Install 'mssearchr' from GitHub:
library(devtools)
install_github("andreysamokhin/mspepsearchr")
```


## Known Discrepancies Between MS Search and MSPepSearch

Users may observe minor differences in results when using MSPepSearch compared
to MS Search. These differences are from the respective NIST tools themselves
and are not caused by this R interface. Any discrepancies can only be addressed
by the developers of the original NIST software.

* Probabilities returned by Identity algorithms may differ slightly (typically
by less than 0.15%).
* Results from the 'Similarity MS/MS in EI' algorithm may vary when the
molecular weight of the unknown compound is automatically estimated from the
mass spectrum.
* When using high-resolution mass spectra, results from algorithms originally
developed for low-resolution spectra (such as 'Identity EI Normal' and
'Similarity EI Simple') may differ.
* RI values are not exported when using the 'Similarity EI Hybrid' algorithm.
* The 'Best matching only' option behaves differently for the 'Identity EI
Normal' and 'Similarity EI Simple' algorithms (MSPepSearch uses only CAS
numbers, whereas MS Search can also use compound names if CAS numbers are
unavailable).


