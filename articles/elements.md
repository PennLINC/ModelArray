# Elements: The Unit of Analysis

## What is an element?

ModelArray performs mass-univariate statistical analysis - fitting the
same model independently at every spatial location in the brain. We call
each spatial location an **element**. Depending on the imaging modality,
an element is one of:

| Modality               | Element type                     | Typical count    | File format            |
|:-----------------------|:---------------------------------|:-----------------|:-----------------------|
| Fixel-based analysis   | Fixel (fiber population element) | ~600,000         | `.mif`                 |
| Voxel-based analysis   | Voxel                            | ~200,000–400,000 | NIfTI (`.nii.gz`)      |
| Surface-based analysis | Greyordinate (cortical vertex)   | ~330,000         | CIFTI (`.dscalar.nii`) |

Despite the differences in what they represent, ModelArray treats all of
these uniformly. The same functions,
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md),
work regardless of whether your elements are fixels, voxels, or
greyordinates.

## The scalar matrix

At the core of every ModelArray object is a **scalar matrix** stored in
an HDF5 file. This matrix has:

- **Rows** = elements (one per spatial location)
- **Columns** = Sources (one per input file)

Each cell contains the scalar value (e.g., FDC, cortical thickness, FA)
for that element in that input file.

``` r
library(ModelArray)

# Create a ModelArray object from an HDF5 file
modelarray <- ModelArray("path/to/data.h5", scalar_types = c("FDC"))

# View the scalar matrix
scalars(modelarray)[["FDC"]]
```

``` console
<602229 x 100> matrix of class DelayedMatrix and type "double":
          FDC/sub-6fee490.mif FDC/sub-647f86c.mif ... FDC/sub-063fd82.mif
     [1,]          0.24264026          0.15679701   .          0.16528498
     [2,]          0.04573315          0.30895054   .          0.33469012
      ...                   .                   .   .                   .
[602229,]          0.12752580          0.51846390   .          0.18338820
```

This is a `DelayedMatrix`, meaning the data lives on disk in the HDF5
file and is only read into memory when accessed. This is what makes
ModelArray memory-efficient — see
[`vignette("hdf5-large-analyses")`](https://pennlinc.github.io/ModelArray/articles/hdf5-large-analyses.md)
for details.

## Element IDs

Every element has an integer **element ID** starting at 0. The element
ID corresponds directly to the row index in the scalar matrix (row 1 =
element ID 0, row 2 = element ID 1, and so on).

When you run a model with ModelArray, the output data frame includes an
`element_id` column that maps results back to their spatial location:

``` r
result <- ModelArray.lm(FDC ~ Age + sex, modelarray, phenotypes, "FDC",
  element.subset = 1:5
)
result$element_id
```

``` console
[1] 0 1 2 3 4
```

You can check the total number of elements with:

``` r
numElementsTotal(modelarray, "FDC")
```

## Working with a subset of elements

The `element.subset` parameter, available in
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
and
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md),
lets you run on a subset of elements rather than all of them. This is
useful for:

- **Testing**: run on the first 100 elements to verify your model works
  before committing to a full run
- **Batch processing**: split a large analysis across multiple HPC jobs
  (see
  [`vignette("element-splitting")`](https://pennlinc.github.io/ModelArray/articles/element-splitting.md))

`element.subset` takes a vector of 1-based indices into the scalar
matrix. For example, `element.subset = 1:100` runs on elements with IDs
0 through 99 (the first 100 rows of the matrix).

``` r
# Test on the first 100 elements
result_test <- ModelArray.lm(FDC ~ Age + sex, modelarray, phenotypes, "FDC",
  element.subset = 1:100
)

# Run on all elements (default)
result_full <- ModelArray.lm(FDC ~ Age + sex, modelarray, phenotypes, "FDC")
```

## Element-specific considerations

For voxel-wise or vertexwise data, different subjects may have different
brain coverage (subject-specific masks). This means some elements near
the edge of the brain may not have valid data from all subjects.

ModelArray handles this with two threshold parameters in
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
and
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md):

- `num.subj.lthr.abs`: minimum number of subjects required (absolute
  count)
- `num.subj.lthr.rel`: minimum proportion of total subjects required
  (0–1)

A voxel is only analyzed if the number of subjects with valid (non-NaN)
values exceeds **both** thresholds. Otherwise, the voxel is skipped and
its outputs are set to `NaN`.

``` r
# Require at least 10 subjects AND at least 20% of total subjects
result <- ModelArray.lm(FA ~ Age + sex, modelarray, phenotypes, "FA",
  num.subj.lthr.abs = 10,
  num.subj.lthr.rel = 0.2
)
```

For fixel-wise and surface-based data, all subjects typically share the
same mask, so these thresholds are not needed.
