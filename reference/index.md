# Package index

## Package overview

- [`ModelArray-package`](https://pennlinc.github.io/ModelArray/reference/ModelArray-package.md)
  : ModelArray: Statistical Analysis of Element-Wise Neuroimaging Data

## ModelArray class and constructor

- [`ModelArray()`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  [`show(`*`<ModelArray>`*`)`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  : ModelArray class
- [`ModelArray()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.md)
  : Load element-wise data from an HDF5 file

## Accessors

- [`sources()`](https://pennlinc.github.io/ModelArray/reference/sources.md)
  : Source filenames of a ModelArray object
- [`scalars()`](https://pennlinc.github.io/ModelArray/reference/scalars.md)
  : Element-wise scalar data of a ModelArray object
- [`results()`](https://pennlinc.github.io/ModelArray/reference/results.md)
  : Statistical results of a ModelArray object
- [`nElements()`](https://pennlinc.github.io/ModelArray/reference/nElements.md)
  : Number of elements in a ModelArray
- [`nInputFiles()`](https://pennlinc.github.io/ModelArray/reference/nInputFiles.md)
  : Number of input files in a ModelArray
- [`scalarNames()`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md)
  : Names of scalars in a ModelArray
- [`analysisNames()`](https://pennlinc.github.io/ModelArray/reference/analysisNames.md)
  : Names of analyses in a ModelArray
- [`elementMetadata()`](https://pennlinc.github.io/ModelArray/reference/elementMetadata.md)
  : Element metadata from a ModelArray
- [`exampleElementData()`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
  : Example per-element data.frame for user functions

## Analysis functions

- [`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
  : Fit element-wise linear models
- [`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
  : Fit element-wise generalized additive models
- [`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
  : Run a user-supplied function for element-wise data

## Per-element internals

These functions are called internally by the analysis functions. Most
users should not call them directly.

- [`analyseOneElement.lm()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.lm.md)
  : Fit a linear model for a single element
- [`analyseOneElement.gam()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.gam.md)
  : Fit a GAM for a single element
- [`analyseOneElement.wrap()`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.wrap.md)
  : Run a user-supplied function for a single element

## Merging

- [`mergeModelArrays()`](https://pennlinc.github.io/ModelArray/reference/mergeModelArrays.md)
  : Merge multiple ModelArrays from different HDF5 files

## Utilities

- [`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)
  : Write outputs from element-wise statistical analysis to an HDF5 file
- [`h5summary()`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
  [`print(`*`<h5summary>`*`)`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
  : Summarize an HDF5 file without loading a full ModelArray
- [`numElementsTotal()`](https://pennlinc.github.io/ModelArray/reference/numElementsTotal.md)
  : Number of elements in ModelArray
- [`gen_gamFormula_fxSmooth()`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_fxSmooth.md)
  : Generate GAM formula with factor-smooth interaction
- [`gen_gamFormula_contIx()`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_contIx.md)
  : Generate GAM formula with continuous-by-continuous interaction
