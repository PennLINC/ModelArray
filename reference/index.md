# Package index

## ModelArray Class

Create, load, and combine ModelArray objects backed by HDF5 files.

- [`ModelArray()`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  [`show(`*`<ModelArray>`*`)`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  : ModelArray class
- [`ModelArray()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.md)
  : Load element-wise data from an HDF5 file
- [`h5summary()`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
  [`print(`*`<h5summary>`*`)`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
  : Summarize an HDF5 file without loading a full ModelArray
- [`mergeModelArrays()`](https://pennlinc.github.io/ModelArray/reference/mergeModelArrays.md)
  : Merge multiple ModelArrays from different HDF5 files

## Accessors

Extract scalar data, source filenames, results, and metadata from a
ModelArray object.

- [`scalars()`](https://pennlinc.github.io/ModelArray/reference/scalars.md)
  : Element-wise scalar data of a ModelArray object
- [`sources()`](https://pennlinc.github.io/ModelArray/reference/sources.md)
  : Source filenames of a ModelArray object
- [`results()`](https://pennlinc.github.io/ModelArray/reference/results.md)
  : Statistical results of a ModelArray object
- [`scalarNames()`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md)
  : Names of scalars in a ModelArray
- [`analysisNames()`](https://pennlinc.github.io/ModelArray/reference/analysisNames.md)
  : Names of analyses in a ModelArray
- [`nElements()`](https://pennlinc.github.io/ModelArray/reference/nElements.md)
  : Number of elements in a ModelArray
- [`nInputFiles()`](https://pennlinc.github.io/ModelArray/reference/nInputFiles.md)
  : Number of input files in a ModelArray
- [`numElementsTotal()`](https://pennlinc.github.io/ModelArray/reference/numElementsTotal.md)
  : Number of elements in ModelArray
- [`elementMetadata()`](https://pennlinc.github.io/ModelArray/reference/elementMetadata.md)
  : Element metadata from a ModelArray

## Element-wise Analysis

Fit statistical models at every element (fixel, voxel, or vertex) in a
ModelArray and write results back to HDF5.

- [`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
  : Fit element-wise linear models

- [`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
  :

  Fit element-wise generalized additive models no model-level p-value
  for GAMs, so there is no `correct.p.value.model` argument.

- [`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
  : Run a user-supplied function for element-wise data

- [`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)
  : Write outputs from element-wise statistical analysis to an HDF5 file

## Testing and Debugging

Build per-element test data and inspect your setup before committing to
a full analysis run.

- [`exampleElementData()`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
  : Example per-element data.frame for user functions

## Formula Helpers

Generate GAM formulas for common interaction structures used with
ModelArray.gam.

- [`gen_gamFormula_fxSmooth()`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_fxSmooth.md)
  : Generate GAM formula with factor-smooth interaction
- [`gen_gamFormula_contIx()`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_contIx.md)
  : Generate GAM formula with continuous-by-continuous interaction

## Internal

Per-element workhorse functions called by the analysis functions. Most
users should not need to call these directly.

- [`analyseOneElement.lm()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.lm.md)
  :

  Fit a linear model for a single element If the number of subjects with
  finite scalar values (not `NaN`, `NA`, or `Inf`) does not exceed
  `num.subj.lthr`, the element is skipped and all statistics are set to
  `NaN`.

- [`analyseOneElement.gam()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.gam.md)
  :

  Fit a GAM for a single element returns metadata (column names, smooth
  term names, parametric term names, and the smoothing parameter
  criterion attribute name) used by `ModelArray.gam` to initialise the
  output data.frame. When `flag_initiate = FALSE`, it returns a numeric
  vector representing one row of the final results matrix.

- [`analyseOneElement.wrap()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.wrap.md)
  : Run a user-supplied function for a single element
