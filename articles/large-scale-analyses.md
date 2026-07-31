# Large-Scale Analyses

## Why HDF5?

A typical fixel-based analysis might involve 600,000 fixels across 1,000
subjects. Stored naively as a double-precision matrix, that’s about 4.5
GB — too large to hold in memory on most systems, let alone alongside
the R session overhead and model fitting.

HDF5 solves this by keeping data on disk and providing efficient random
access to slices of the data. ModelArray leverages this to fit models
element-by-element without ever loading the full matrix into RAM.

**TileDB**

TileDB solves the same on-disk storage problem with a directory-backed
`.tdb` store instead of a single `.h5` file. ModelArray uses the same
delayed access pattern for TileDB-backed stores, so the modelling
workflow stays the same once the store is loaded.

## Chunking

HDF5 files store data in **chunks** — rectangular blocks of the array
that are read and written as a unit. When you access a single row
(element), HDF5 only reads the chunk(s) that contain that row, not the
entire dataset.

The chunk layout matters for performance. Because ModelArray reads data
**row by row** (one element at a time), chunks that span a modest number
of rows and all columns work best. This is the default behavior of the
`modelarrayio to-modelarray` conversion command.

When creating HDF5 files with `modelarrayio to-modelarray`, you can
control the target chunk size with the `--target-chunk-mb` flag:

``` console
$ modelarrayio to-modelarray \
    --backend hdf5 \
    --index-file FDC/index.mif \
    --directions-file FDC/directions.mif \
    --cohort-file cohort.csv \
    --output data.h5 \
    --target-chunk-mb 4       # target chunk size in MiB
```

**TileDB**

Use `--backend tiledb` and a `.tdb` output path. In TileDB, the
storage-level analogue of an HDF5 chunk is a **tile**. ModelArrayIO uses
the same `--target-chunk-mb` option to choose an approximate TileDB tile
size:

``` console
$ modelarrayio to-modelarray \
    --backend tiledb \
    --index-file FDC/index.mif \
    --directions-file FDC/directions.mif \
    --cohort-file cohort.csv \
    --output data.tdb \
    --target-chunk-mb 4       # target TileDB tile size in MiB
```

This controls the on-disk TileDB layout. It is separate from
ModelArray’s analysis-time TileDB read-block settings, described below
under memory efficiency.

For HDF5, smaller chunks use less memory per read but may require more
disk seeks. The default target is 2 MiB; the examples above explicitly
override it to 4 MiB.

## Compression

HDF5 supports transparent compression of chunks using gzip. When
enabled, data is compressed on disk and decompressed on-the-fly when
read. This typically reduces file sizes by 30–60% with minimal impact on
read speed.

`ModelArrayIO` enables gzip compression by default. You can see this
reflected in the output filenames:

``` console
study-HBN_compression-gzip_chunkmb-4_thickness.h5
```

The compression is transparent to ModelArray — it reads the data the
same way regardless of whether it’s compressed.

**TileDB**

TileDB stores are directories, so the output path ends in `.tdb`:

``` console
study-HBN_compression-gzip_chunkmb-4_thickness.tdb
```

ModelArrayIO’s compression options map to TileDB compression filters and
tile sizes. As with HDF5, compression is transparent when you load the
store in R.

## Memory efficiency during model fitting

When you call
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
or
[`ModelArray.wrap()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md),
here is what happens for each element:

1.  **Read one row** from the scalar matrix (the values for all subjects
    at this element)
2.  **Combine** with the phenotypes data frame to create a per-element
    data frame
3.  **Fit the model** (lm, gam, or your custom function)
4.  **Extract statistics** into a single row of the output data frame
5.  **Discard** the per-element data and move to the next element

At no point is more than one row of the scalar matrix in memory. This is
why ModelArray can handle arbitrarily large datasets — the memory
footprint is determined by the number of subjects (columns), not the
number of elements (rows).

**TileDB**

For TileDB-backed data, ModelArray may hold a bounded block of rows in
memory to avoid slow one-row TileDB reads. The memory footprint is still
controlled: it depends on the number of subjects, the number of attached
scalars, and the configured block size, not the full number of elements
in the store.

The default TileDB read-block target is 512 MiB. Tune it when needed:

``` r

# Target roughly 256 MiB per TileDB read block
options(ModelArray.tiledb_read_block_mb = 256)

# Or set an explicit number of elements per block
options(ModelArray.tiledb_read_block_size = 5000)

# Disable TileDB row-block materialization and read row by row
options(ModelArray.tiledb_read_block_size = 0)
```

Approximate memory for one TileDB read block is:

``` text
block_size x number_of_subjects x number_of_attached_scalars x 8 bytes
```

For example, a block of 5,000 elements across 1,000 subjects for one
scalar is about 40 MB before R object overhead.

## Parallelism within a single process

ModelArray supports parallel processing via the `n_cores` parameter:

``` r

result <- ModelArray.lm(FDC ~ Age + sex, modelarray, phenotypes, "FDC",
  n_cores = 4
)
```

This uses
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(fork-based parallelism on Linux/macOS) to process multiple elements
simultaneously. Each forked worker inherits the same file handle and
reads different rows from the HDF5 file. Because HDF5 supports
concurrent **reads** from the same file, this works safely and scales
well.

**TileDB**

For TileDB-backed scalars, ModelArray reads each row block in the parent
process before dispatching elements to workers. On Linux/macOS this lets
forked workers reuse the in-memory block through copy-on-write, reducing
repeated TileDB reads.

## Why concurrent writes don’t work

While multiple processes can safely **read** from the same HDF5 file,
they **cannot safely write** to it simultaneously. The standard HDF5
library (without special MPI-IO compilation) does not support concurrent
writes — doing so can corrupt the file.

This means you should **never** have multiple R processes calling
[`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)
on the same HDF5 file at the same time. If you split your analysis
across HPC jobs (see
[`vignette("element-splitting")`](https://pennlinc.github.io/ModelArray/articles/element-splitting.md)),
each job should save its partial results to a separate `.rds` file, then
a single final step combines them and writes to the HDF5 file.

``` r

# WRONG: multiple jobs writing to the same H5 simultaneously
# writeResults("data.h5", df.output = my_partial_result, ...)  # Don't do this in parallel!

# RIGHT: save partial results, then combine in one process
saveRDS(my_partial_result, sprintf("result_chunk_%d.rds", job_id))

# ... after all jobs finish, in a single process:
chunks <- lapply(Sys.glob("result_chunk_*.rds"), readRDS)
full_result <- do.call(rbind, chunks)
writeResults("data.h5", df.output = full_result, analysis_name = "results_lm")
```

**TileDB**

Use the same single-writer pattern for ModelArray writes to TileDB
stores. Even though TileDB has its own concurrency features,
[`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)
is intended to create or replace one analysis result from one R process
at a time.

``` r

# WRONG: multiple jobs writing to the same TileDB store simultaneously
# writeResults("data.tdb", df.output = my_partial_result, backend = "tiledb")

# RIGHT: save partial results, then combine in one process
saveRDS(my_partial_result, sprintf("result_chunk_%d.rds", job_id))

# ... after all jobs finish, in a single process:
chunks <- lapply(Sys.glob("result_chunk_*.rds"), readRDS)
full_result <- do.call(rbind, chunks)
writeResults("data.tdb",
  df.output = full_result,
  analysis_name = "results_lm",
  backend = "tiledb"
)
```

## Backing up your HDF5 file

Because
[`writeResults()`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)
modifies the HDF5 file in place, it’s good practice to back up your file
before writing results:

``` console
$ cp data.h5 data_backup.h5
```

If something goes wrong during a write (e.g., the process is killed
mid-write), the HDF5 file may be left in an inconsistent state. Having a
backup lets you recover without re-running the conversion step.

**TileDB**

TileDB stores are directories, so back up the full directory before
writing results:

``` console
$ cp -R data.tdb data_backup.tdb
```

If a TileDB write is interrupted, restore from the backup directory
before rerunning the write step.
