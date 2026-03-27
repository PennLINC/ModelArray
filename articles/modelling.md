# Modelling in Depth: lm, gam, and wrap

This vignette covers the three modelling functions in ModelArray in
detail. If you haven’t already, start with
[`vignette("walkthrough")`](https://pennlinc.github.io/ModelArray/articles/walkthrough.md)
for the end-to-end workflow, and
[`vignette("elements")`](https://pennlinc.github.io/ModelArray/articles/elements.md)
for background on what an element is.

``` r
library(ModelArray)

# Using the PNC fixel demo data (see walkthrough for download instructions)
h5_path <- "~/Desktop/myProject/demo_FDC_n100.h5"
modelarray <- ModelArray(h5_path, scalar_types = c("FDC"))
phenotypes <- read.csv("~/Desktop/myProject/cohort_FDC_n100.csv")
```

## Linear models with `ModelArray.lm()`

### Basic usage

[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
fits [`lm()`](https://rdrr.io/r/stats/lm.html) independently at every
element. Define a formula with the scalar name as the response variable:

``` r
formula.lm <- FDC ~ Age + sex + dti64MeanRelRMS

result <- ModelArray.lm(formula.lm, modelarray, phenotypes, "FDC",
  element.subset = 1:100
)
```

### Output columns

By default, for each term (Intercept, Age, sex, etc.) the output
includes:

| Column               | Description                  |
|:---------------------|:-----------------------------|
| `<term>.estimate`    | Coefficient estimate (slope) |
| `<term>.statistic`   | *t*-statistic                |
| `<term>.p.value`     | Uncorrected *p*-value        |
| `<term>.p.value.fdr` | FDR-corrected *p*-value      |

And for the overall model:

| Column                | Description                   |
|:----------------------|:------------------------------|
| `model.adj.r.squared` | Adjusted *R*-squared          |
| `model.p.value`       | Model *F*-test *p*-value      |
| `model.p.value.fdr`   | FDR-corrected model *p*-value |

### More comprehensive output

Use `full.outputs = TRUE` to get additional statistics:

``` r
result_full <- ModelArray.lm(formula.lm, modelarray, phenotypes, "FDC",
  element.subset = 1:100,
  full.outputs = TRUE
)
```

### P-value correction

By default, FDR correction is applied. You can request additional
correction methods:

``` r
result <- ModelArray.lm(formula.lm, modelarray, phenotypes, "FDC",
  element.subset = 1:100,
  correct.p.value.terms = c("fdr", "bonferroni")
)
```

This adds `<term>.p.value.bonferroni` columns alongside the
FDR-corrected ones.

## Generalized Additive Models with `ModelArray.gam()`

GAMs are useful when you expect nonlinear relationships - for example,
nonlinear changes in brain structure across the lifespan.

### Basic usage

Use `s()` to specify smooth terms. The `k` parameter sets an upper limit
on the degrees of freedom for the smooth function:

``` r
formula.gam <- FDC ~ s(Age, k = 4) + sex + dti64MeanRelRMS

result <- ModelArray.gam(formula.gam, modelarray, phenotypes, "FDC",
  element.subset = 1:100,
  method = "REML"
)
```

ModelArray prints a summary of the smooth term settings when you run
this:

``` console
The formula requested: FDC ~ s(Age, k = 4) + sex + dti64MeanRelRMS
s(Age):   k = 4;   fx = FALSE (default);   bs = tp (default)
method = REML (default: GCV.Cp)
```

### Output columns

GAM output differs from lm for smooth terms vs parametric terms:

**Smooth terms** (e.g., `s(Age)`):

| Column              | Description                  |
|:--------------------|:-----------------------------|
| `s_Age.statistic`   | *F*-statistic for the smooth |
| `s_Age.p.value`     | Uncorrected *p*-value        |
| `s_Age.p.value.fdr` | FDR-corrected *p*-value      |

**Parametric terms** (e.g., `sex`, `dti64MeanRelRMS`):

| Column             | Description           |
|:-------------------|:----------------------|
| `<term>.estimate`  | Coefficient estimate  |
| `<term>.statistic` | *t*-statistic         |
| `<term>.p.value`   | Uncorrected *p*-value |

**Model-level**:

| Column           | Description                           |
|:-----------------|:--------------------------------------|
| `model.dev.expl` | Proportion of null deviance explained |

### P-value correction for GAMs

Smooth terms and parametric terms have separate correction arguments:

``` r
result <- ModelArray.gam(formula.gam, modelarray, phenotypes, "FDC",
  element.subset = 1:100,
  correct.p.value.smoothTerms = c("fdr", "bonferroni"),
  correct.p.value.parametricTerms = c("fdr", "bonferroni"),
  method = "REML"
)
```

### Smooth term options

ModelArray supports the full range of
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) smooth terms:

- `s()`: penalized regression splines (default thin-plate, `bs = "tp"`)
- `te()`, `ti()`, `t2()`: tensor product smooths for interactions
  between continuous variables

The `fx` parameter controls whether the smooth is penalized
(`fx = FALSE`, the default) or unpenalized/fixed (`fx = TRUE`). The `bs`
parameter selects the basis type.

### GAM formula helpers

ModelArray provides two helper functions for common GAM formula
patterns:

**Factor-smooth interaction** (`gen_gamFormula_fxSmooth`): Generates
`y ~ orderedFactor + s(x) + s(x, by=orderedFactor)`:

``` r
result <- gen_gamFormula_fxSmooth(
  response.var = "FDC",
  factor.var = "group",
  smooth.var = "Age",
  phenotypes = phenotypes,
  reference.group = "control",
  fx = TRUE, k = 4
)
formula <- result$formula
phenotypes <- result$phenotypes # may contain new ordered factor column
```

**Continuous interaction** (`gen_gamFormula_contIx`): Generates
`y ~ ti(x) + ti(z) + ti(x, z)`:

``` r
formula <- gen_gamFormula_contIx(
  response.var = "FDC",
  cont1.var = "Age",
  cont2.var = "BMI",
  fx = TRUE, k = 4
)
```

## Custom functions with `ModelArray.wrap()`

[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
lets you run any R function across all elements, using ModelArray’s
efficient looping and parallelization machinery.

### Prototyping with `exampleElementData()`

Before writing your function, use
[`exampleElementData()`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
to see what a per-element data frame looks like:

``` r
# Using the HBN cortical thickness data
h5_path <- "path/to/thickness.h5"
modelarray <- ModelArray(h5_path, scalar_types = c("thickness"))
phenotypes <- read.csv("path/to/phenotypes.csv")

example_df <- exampleElementData(modelarray,
  scalar = "thickness", i_element = 12345, phenotypes = phenotypes
)
names(example_df)
```

This returns a data frame with all the phenotype columns plus a new
column for the scalar value at that element. You can explore it like any
data frame:

``` r
library(ggplot2)
ggplot(example_df, aes(x = age, y = thickness, color = study_site)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_minimal()
```

### Writing your function

Your function must accept a single argument named `data` and return a
single-row data frame (or tibble):

``` r
my_site_analysis <- function(data) {
  # Compute mean thickness per site
  site_means <- data %>%
    dplyr::group_by(study_site) %>%
    dplyr::summarize(mean_thickness = mean(thickness, na.rm = TRUE)) %>%
    tidyr::pivot_wider(
      names_from = study_site,
      values_from = mean_thickness,
      names_glue = "thickness_{study_site}"
    )

  # Run ANOVA on site
  site_anova <- lm(thickness ~ study_site, data = data) %>%
    anova() %>%
    broom::tidy() %>%
    dplyr::filter(term != "Residuals") %>%
    dplyr::select(term, statistic, p.value) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(statistic, p.value),
      names_glue = "{term}.{.value}"
    )

  dplyr::bind_cols(site_means, site_anova)
}
```

### Testing on a single element

Use
[`analyseOneElement.wrap()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.wrap.md)
to verify your function works in ModelArray’s machinery:

``` r
analyseOneElement.wrap(
  12345, my_site_analysis, modelarray, phenotypes, "thickness",
  num.subj.lthr = 10, flag_initiate = TRUE, on_error = "debug"
)
```

### Running across all elements

``` r
result <- ModelArray.wrap(
  FUN = my_site_analysis,
  data = modelarray,
  phenotypes = phenotypes,
  scalar = "thickness",
  n_cores = 8
)

head(result)
writeResults(h5_path, df.output = result, analysis_name = "site_analysis")
```

### Case study: harmonize with `covfam` and stream to a new scalar

You can also use
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
as a transformation engine rather than a test runner. The example below
applies a harmonization step per element, streams the harmonized values
directly into `/scalars`, and then loads that harmonized scalar in a new
`ModelArray`.

``` r
library(covfam)

# Assume phenotypes has harmonization variables, e.g.:
#   source_file, site, age, sex
# and that "thickness" is already present in the h5 as an input scalar.

h5_in <- "thickness_raw.h5"
h5_harmonized <- "thickness_harmonized.h5"

ma_raw <- ModelArray(h5_in, scalar_types = "thickness")
phenotypes <- read.csv("phenotypes.csv")

# Copy metadata/layout into a new file so we can append new scalars there.
file.copy(h5_in, h5_harmonized, overwrite = TRUE)

covfam_harmonize_element <- function(data) {
  # Harmonize one element across subjects.
  # Adjust arguments here to match your covfam configuration.
  out <- covfam::covfam(
    y = data$thickness,
    batch = data$site,
    mod = data.frame(age = data$age, sex = data$sex)
  )

  # Return named subject-level vector so wrap columns map to source_file order.
  vals <- out$y_harmonized
  names(vals) <- data$source_file
  vals
}

# Stream harmonized subject-level values into scalars/thickness_covfam/values.
# Set return_output = FALSE to avoid keeping the full output table in memory.
ModelArray.wrap(
  FUN = covfam_harmonize_element,
  data = ma_raw,
  phenotypes = phenotypes,
  scalar = "thickness",
  n_cores = 8,
  write_scalar_name = "thickness_covfam",
  write_scalar_file = h5_harmonized,
  write_scalar_flush_every = 2000L,
  return_output = FALSE
)

# Load the harmonized scalar as a new modelling input
ma_harmonized <- ModelArray(h5_harmonized, scalar_types = c("thickness_covfam"))

# Use harmonized data in downstream models
fit_harmonized <- ModelArray.lm(
  thickness_covfam ~ age + sex,
  data = ma_harmonized,
  phenotypes = phenotypes,
  scalar = "thickness_covfam",
  n_cores = 8
)
```

If you want to stream model statistics too (not only transformed
scalars), use `write_results_name` and `write_results_file` in
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
or
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md).

#### Quick QA: export one harmonized image to inspect

After writing `thickness_covfam` into `/scalars`, you can export one row
back to an image and open it in your usual viewer as a spot-check. Use
`--column-index` to select which row from
`scalars/thickness_covfam/values` to export.

``` bash
modelarrayio h5-export-nifti-file \
  --input-hdf5 thickness_harmonized.h5 \
  --scalar-name thickness_covfam \
  --column-index 0 \
  --group-mask-file group_mask_thickness.nii.gz \
  --output-file qa_thickness_covfam_col000.nii.gz
```

If your workflow uses CIFTI or MIF instead of NIfTI, use the matching
export command: `h5-export-cifti-file` or `h5-export-mif-file`.

## Modelling across multiple h5 files with `mergeModelArrays()`

When scalars live in separate h5 files — for example, cortical thickness
in one file and curvature in another — you can combine them for
cross-scalar modelling without rewriting or merging h5 files on disk.

### Setup

Each h5 file gets its own ModelArray and phenotypes data frame. The
phenotypes must share a common identifier column (e.g., `subject_id`):

``` r
ma_thick <- ModelArray("thickness.h5", scalar_types = "thickness")
ma_curv <- ModelArray("curv.h5", scalar_types = "curv")

phen_thick <- read.csv("phenotypes_thickness.csv")
# Must contain: subject_id, source_file, (other covariates)

phen_curv <- read.csv("phenotypes_curv.csv")
# Must contain: subject_id, source_file
```

### Merging

``` r
merged <- mergeModelArrays(
  list(ma_thick, ma_curv),
  list(phen_thick, phen_curv),
  merge_on = "subject_id"
)
```

This performs an inner join on `subject_id`, keeping only subjects
present in both files. The result is a list with:

- `merged$data` — a combined ModelArray with both `thickness` and `curv`
  scalars
- `merged$phenotypes` — the joined phenotypes, with
  `source_file.thickness` and `source_file.curv` columns preserving the
  original filenames, plus a unified `source_file` column used
  internally

### Cross-scalar modelling

Now you can use one scalar as a predictor for another:

``` r
result <- ModelArray.lm(
  thickness ~ curv + Age + sex,
  data = merged$data,
  phenotypes = merged$phenotypes,
  scalar = "thickness"
)
```

This works with
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
and
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
too.

### N-way merges

[`mergeModelArrays()`](https://pennlinc.github.io/ModelArray/reference/mergeModelArrays.md)
accepts any number of ModelArrays:

``` r
merged <- mergeModelArrays(
  list(ma_thick, ma_curv, ma_sulc),
  list(phen_thick, phen_curv, phen_sulc),
  merge_on = "subject_id"
)
```

All scalar names must be unique across the inputs.

## Converting results back to image format

After
[`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md),
use the appropriate ModelArrayIO tool to convert results back to the
original image format for visualization:

| Data type | Command            | Viewer               |
|:----------|:-------------------|:---------------------|
| Fixel     | `fixelstats_write` | MRtrix MRView        |
| Voxel     | `voxelstats_write` | FSLeyes, MRView      |
| Surface   | `ciftistats_write` | Connectome Workbench |

See the [ModelArrayIO
documentation](https://github.com/PennLINC/ModelArrayIO) for
command-line usage details, and
[`vignette("walkthrough")`](https://pennlinc.github.io/ModelArray/articles/walkthrough.md)
for a worked example with MRView.
