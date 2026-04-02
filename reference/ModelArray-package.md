# ModelArray: Statistical Analysis of Element-Wise Neuroimaging Data

The ModelArray package provides an S4 class and associated methods for
performing massively univariate statistical analyses on element-wise
(fixel, voxel, or vertex) neuroimaging data stored in HDF5 files.

## Details

The core workflow is:

1.  Inspect an HDF5 file with
    [`h5summary()`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)

2.  Load data with
    [`ModelArray()`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)

3.  Fit models with
    [`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
    [`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
    or
    [`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)

4.  Access results via
    [`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
    [`results`](https://pennlinc.github.io/ModelArray/reference/results.md),
    [`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md)

For multi-session or multi-modality analyses, combine ModelArrays with
[`mergeModelArrays`](https://pennlinc.github.io/ModelArray/reference/mergeModelArrays.md)
before fitting cross-scalar models.

## Centralized imports

The following imports are consolidated here because they cannot be
expressed as `pkg::fun()` qualified calls. All other add-on package
functions are called with explicit namespacing throughout the package
source code.

## References

Zhao, C., et al. (2023). ModelArray–An R package for statistical
analysis of fixel-wise data. *NeuroImage*, 271, 120037.
[doi:10.1016/j.neuroimage.2023.120037](https://doi.org/10.1016/j.neuroimage.2023.120037)

## See also

[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md),
[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html),
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md),
[`mergeModelArrays`](https://pennlinc.github.io/ModelArray/reference/mergeModelArrays.md),
[`h5summary`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)

## Author

**Maintainer**: Matthew Cieslak <Matthew.Cieslak@pennmedicine.upenn.edu>

Authors:

- Chenying Zhao

- Tinashe Tapera

- Theodore Satterthwaite
