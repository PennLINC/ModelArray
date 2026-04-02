# Merge multiple ModelArrays from different HDF5 files

Combines scalars from multiple
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
objects into a single
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md),
aligning subjects via shared phenotype columns.

## Usage

``` r
mergeModelArrays(modelarrays, phenotypes_list, merge_on)
```

## Arguments

- modelarrays:

  A list of at least two
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  objects, each constructed from a different HDF5 file.

- phenotypes_list:

  A list of data.frames, one per
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  in `modelarrays`. Each must contain a `source_file` column whose
  entries match the corresponding ModelArray's sources (i.e.
  `sources(modelarrays[[i]])`). Each must also contain all columns named
  in `merge_on`.

- merge_on:

  Character vector of column names present in all data.frames in
  `phenotypes_list`, used to inner-join subjects across
  sessions/modalities (e.g. `c("subject_id")`). The combination of these
  columns must uniquely identify each subject within each data.frame.

## Value

A list with two components:

- data:

  A combined
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  containing scalars from all inputs. Each scalar's columns are
  subsetted and reordered to match the inner-joined subject list.

- phenotypes:

  The inner-joined data.frame. Original `source_file` columns are
  renamed to `source_file.<scalar_name>` and a new unified `source_file`
  column is added for use with analysis functions.

## Details

The merge performs an inner join of the phenotype data.frames on the
columns specified by `merge_on`. Only subjects present in all phenotype
data.frames are retained. Scalar matrices from each input
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
are column-subsetted and reordered to match the joined subject list.

A unified `source_file` column is created from the `merge_on` columns so
that downstream analysis functions
([`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md))
can align phenotypes to scalars. The original `source_file` columns are
renamed to `source_file.<first_scalar_name>` for each input
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).

Scalar names must be unique across all input ModelArrays. If two
ModelArrays share a scalar name (e.g. both have `"FD"`), the function
will error. Element counts (number of rows) must match across all
scalars.

If element metadata is available (see
[`elementMetadata`](https://pennlinc.github.io/ModelArray/reference/elementMetadata.md)),
the function checks that it is consistent across inputs and warns if it
differs or is only partially available.

## See also

[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for constructing individual ModelArray objects,
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
for fitting models on the merged object,
[`elementMetadata`](https://pennlinc.github.io/ModelArray/reference/elementMetadata.md)
for element correspondence checks.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load two sessions from different h5 files
ma1 <- ModelArray("session1.h5", scalar_types = c("FD"))
ma2 <- ModelArray("session2.h5", scalar_types = c("FC"))
phen1 <- read.csv("session1_cohort.csv")
phen2 <- read.csv("session2_cohort.csv")

# Merge on subject ID
merged <- mergeModelArrays(
  modelarrays = list(ma1, ma2),
  phenotypes_list = list(phen1, phen2),
  merge_on = "subject_id"
)

# Use the merged object for cross-scalar analysis
merged$data
scalarNames(merged$data) # c("FD", "FC")
head(merged$phenotypes)

results <- ModelArray.lm(
  FD ~ age + sex + FC,
  data = merged$data,
  phenotypes = merged$phenotypes,
  scalar = c("FD", "FC")
)
} # }
```
