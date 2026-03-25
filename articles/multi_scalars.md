# Using Multiple Scalars in ModelArray

``` r
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

### Overview

ModelArray supports data sets that contain more than one element-wise
scalar (e.g., FA and MD). This vignette shows how to:

- Fit models where the response is one scalar and predictors can include
  other scalars.
- Ensure sources alignment and avoid name collisions in `phenotypes`.
- Use `ModelArray.wrap` to run custom per-element functions with all
  scalars attached.

We will use both the package example data and small in-memory examples
for clarity and speed.

### Setup

``` r
library(ModelArray)
library(dplyr)
set.seed(123)

h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
phenotypes <- read.csv(csv_path)
```

### In-memory multi-scalar object

``` r
# Build a tiny in-memory ModelArray with two scalars (FA, MD)
num_elements <- 3L
num_subj <- nrow(phenotypes)
subj <- phenotypes$source_file

MD <- matrix(rnorm(num_elements * num_subj), nrow = num_elements, ncol = num_subj)
FA <- matrix(NA_real_, nrow = num_elements, ncol = num_subj)
for (i in seq_len(num_elements)) {
  FA[i, ] <- 0.6 * MD[i, ] + 0.05 * phenotypes$age + rnorm(num_subj, sd = 0.05)
}

ma <- methods::new(
  "ModelArray",
  scalars = list(FA = FA, MD = MD),
  sources = list(FA = subj, MD = subj),
  results = list(),
  path = ""
)
```

### Linear models with scalar predictors

You can reference another scalar (here `MD`) as a predictor when
modeling `FA`:

``` r
res_lm <- ModelArray.lm(
  FA ~ MD + age,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FA",
  element.subset = as.integer(1:2),
  num.subj.lthr.abs = 2L, num.subj.lthr.rel = 0,
  verbose = FALSE, pbar = FALSE
)
colnames(res_lm)[1:10]
```

Requirements and helpful checks:

- The formula response (LHS) must match `scalar` (here `FA`).
- `phenotypes$source_file` is aligned to `sources(data)[[scalar]]`
  automatically, with informative errors if mismatched.
- If `phenotypes` contains a column named like any scalar (e.g., `MD`),
  modeling will error with a clear collision message. Rename such
  columns.

### GAMs with scalar predictors

Similarly, you may use scalar predictors in smooth terms:

``` r
res_gam <- ModelArray.gam(
  FA ~ s(MD) + s(age),
  data = ma,
  phenotypes = phenotypes,
  scalar = "FA",
  element.subset = as.integer(1:2),
  num.subj.lthr.abs = 2L, num.subj.lthr.rel = 0,
  verbose = FALSE, pbar = FALSE
)
colnames(res_gam)[1:10]
```

The same alignment and collision rules apply. If the predictor scalar’s
source list does not match `phenotypes$source_file`, you will get an
informative error describing the mismatch.

### Wrapping custom functions with all scalars

`ModelArray.wrap` attaches all scalars into the per-element data.frame
before calling your function. This enables custom analyses that use
multiple scalars at once.

``` r
my_fun <- function(data) {
  stopifnot(all(c("FA", "MD") %in% names(data)))
  # Example: simple correlation between FA and MD for the current element
  keep <- is.finite(data$FA) & is.finite(data$MD)
  r <- suppressWarnings(stats::cor(data$FA[keep], data$MD[keep]))
  tibble::tibble(cor_FA_MD = r)
}

res_wrap <- ModelArray.wrap(
  FUN = my_fun,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FA",
  element.subset = as.integer(1:2),
  num.subj.lthr.abs = 2L, num.subj.lthr.rel = 0,
  verbose = FALSE, pbar = FALSE
)
res_wrap
```

Notes:

- `ModelArray.wrap` checks for collisions between scalar names and
  `phenotypes` columns and errors with guidance to rename.
- Sources must match for every scalar attached; mismatches generate
  informative errors identifying the offending scalar.
- Subject-thresholding uses the intersection of valid (finite) values
  across all attached scalars.

### Tips for avoiding common pitfalls

- If you need a `phenotypes` column with the same name as a scalar
  (e.g., for an external measurement), rename it (e.g., `MD_ext`).
- When bringing additional scalars, ensure their `sources` vectors refer
  to the same subject IDs as `phenotypes$source_file` (order is handled
  automatically).

### See also

- The wrapping vignette: see `vignettes/wrap_function.Rmd` for a broader
  introduction to custom per-element functions.

``` r
# For pkgdown sites, this will appear in the Articles menu
```
