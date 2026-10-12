# R analysis

Run after peak calling. Add your R scripts here, for example:

- `01_permutation_tests.R` - regioneR permutation tests (numOverlaps for
  spKAS-seq; meanDistance for the DIvA RIAN-seq libraries), with
  Benjamini-Hochberg correction across the testable samples.
- `02_annotation.R` - ChIPseeker annotation against the CHM13v2.0 GTF, with
  chrM peaks removed and promoters defined as within 3 kb of a TSS.
- `03_figures.R` - figure code.

Save them from RStudio before uploading. If some analysis was run in the
console, use the History tab ("To Source") to recover it.
