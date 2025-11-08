# mspepsearchr

## mspepsearchr 0.2.0

* Added external process parallelization using the `parallel` package.
* Created the "Introduction to mspepsearchr" vignette.
* Changed the default value of several arguments:
  * `load_in_memory` - `TRUE`;
  * `precursor_ion_tol` - `list(value = 20, utits = "ppm")`;
  * `product_ions_tol` - `list(value = 0.01, units = "mz")`;
  * `ignore_precursor_ion_tol` - `list(value = 1.6, units = "mz")`.
* Added the *fenamiphos_msms_hr.msp* file.
* Updated the *massbank_subset_msms_hr* library (added three compounds with
precursor m/z 304).
* Refactored internal functions.
* Updated and expanded test coverage.


## mspepsearchr 0.1.0

* The initial release.
