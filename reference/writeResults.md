# Write outputs from element-wise statistical analysis to the HDF5 file.

Create a group named \`analysis_name\` in HDF5 file, then write the
statistical results data.frame (i.e. for one analysis) in it.

## Usage

``` r
writeResults(
  fn.output,
  df.output,
  analysis_name = "myAnalysis",
  overwrite = TRUE
)
```

## Arguments

- fn.output:

  A character, The HDF5 (.h5) filename for the output

- df.output:

  A data.frame object with element-wise statistical results, returned
  from \`ModelArray.lm()\` etc

- analysis_name:

  A character, the name of the results

- overwrite:

  If a group with the same analysis_name exists in HDF5 file, whether
  overwrite it (TRUE) or not (FALSE)

## Details

debug tip: For "Error in H5File.open(filename, mode, file_create_pl,
file_access_pl)", check if there is message 'No such file or directory'.
Try absolute .h5 filename.
