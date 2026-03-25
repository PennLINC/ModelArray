# Application to voxel-wise data

`ModelArray` is compatible with not only fixel-wise data, but also
voxel-wise data. The
[`vignette("walkthrough")`](https://pennlinc.github.io/ModelArray/articles/walkthrough.md)
page has shown example walkthrough with example fixel-wise data, with
hints for voxel-wise data provided throughout the walkthrough. In
general, the workflow of applying `ModelArray` to voxel-wise data is
consistent to that of fixel-wise data, but please note several important
differences when execution described as below.

## Conversion of voxel-wise data format

When converting voxel-wise data between NIfTI and HDF5 file format that
`ModelArray` requires, you will still use `ConFixel` software, but a
different converter tailored for voxel-wise data, `ConVoxel`. We have
provided a walkthrough of how to use it
[here](https://github.com/PennLINC/ConFixel/blob/main/notebooks/walkthrough_voxel-wise_data.md).

## Applying `ModelArray` to voxel-wise data

You can apply the same `ModelArray` functions, arguments, etc used in
fixel-wise data to voxel-wise data; just be aware that an “element” is a
“voxel” now.

In addition, there are additional arguments specifically useful for
voxel-wise data, regarding the **subject-specific masks** (i.e.,
individual masks). Different from fixel-wise data where subjects usually
share the same mask, voxel-wise data may have *subject-specific masks*,
where different subjects may have different brain coverage even in
template space. Therefore, for voxels (e.g., at the edge of the brain)
which do not have data from sufficient number of subjects, their
statistic results might not be reliable. It’s better to exclude those
voxels from the analysis.

In `ConVoxel`, for voxels not within a subject-specific mask (but within
a larger, group-level mask), their values will be set to `NaN`. Thus
`ModelArray` can know the number of subjects that have valid values in
each voxel within the group mask.

To set up minimal number of subjects required for a voxel, you can use
these two arguments in
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
and
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md):

- `num.subj.lthr.abs`: lower threshold, absolute number of subjects;
  e.g., 10
- `num.subj.lthr.rel`: lower threshold, proportion of total number of
  subjects, a value between 0-1; e.g., 0.2 (i.e., 20% of all subjects)

Suppose a voxel has valid values from N subjects, where N:

- N \> `num.subj.lthr.abs`, and
- N \> `num.subj.lthr.rel` \* `Ntotal`

where `Ntotal` is the total number of subjects, which is quantified by
number of rows in argument `phenotypes`. If N satisfied above criteria,
statistical analysis of this voxel will be performed as normal;
otherwise, this voxel will be skipped and statistical outputs will be
set as `NaN`.

For more details about these two arguments, please check the reference
pages of
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
and
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md).

Finally, if you choose to save number of observations used (i.e., if
specifying `nobs` in `var.model` in
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
and
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)),
after conversion using `ConVoxel`, you’ll see an image called
`*_model.nobs*`. If the subject-specific masks you provided are
different across subjects, in this image you’ll probably see the
heterogeneity of number of observations used in the brain, with fewer
observations used at the edge of the brain.
