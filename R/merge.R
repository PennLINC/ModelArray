#' Merge multiple ModelArrays from different h5 files
#'
#' @description
#' Combines scalars from multiple ModelArray objects into a single ModelArray,
#' aligning subjects via phenotype columns. Uses `DelayedArray::acbind()` for
#' virtual column-binding — no h5 rewriting is needed.
#'
#' @param modelarrays A list of ModelArray objects
#' @param phenotypes_list A list of data.frames, one per ModelArray. Each must
#'   contain a `source_file` column matching its corresponding ModelArray's sources.
#' @param merge_on Character vector of column names to join phenotypes on
#'   (e.g., `c("subject_id")`). These columns must exist in all phenotypes data.frames.
#' @return A list with:
#'   \item{data}{A combined ModelArray with scalars from all inputs}
#'   \item{phenotypes}{The inner-joined phenotypes data.frame. The original
#'     `source_file` columns are renamed to `source_file.<scalar_name>`.}
#' @importFrom DelayedArray acbind
#' @importFrom utils head
#' @export
mergeModelArrays <- function(modelarrays, phenotypes_list, merge_on) {
  if (!is.list(modelarrays) || length(modelarrays) < 2) {
    stop("modelarrays must be a list of at least 2 ModelArray objects")
  }
  if (length(modelarrays) != length(phenotypes_list)) {
    stop("modelarrays and phenotypes_list must have the same length")
  }

  # Validate inputs
  for (i in seq_along(modelarrays)) {
    if (!inherits(modelarrays[[i]], "ModelArray")) {
      stop("modelarrays[[", i, "]] is not a ModelArray object")
    }
    if (!is.data.frame(phenotypes_list[[i]])) {
      stop("phenotypes_list[[", i, "]] is not a data.frame")
    }
    if (!("source_file" %in% colnames(phenotypes_list[[i]]))) {
      stop("phenotypes_list[[", i, "]] must contain a 'source_file' column")
    }
    for (col in merge_on) {
      if (!(col %in% colnames(phenotypes_list[[i]]))) {
        stop("Column '", col, "' not found in phenotypes_list[[", i, "]]")
      }
    }
  }

  # Check for scalar name collisions
  all_scalar_names <- unlist(lapply(modelarrays, function(ma) names(scalars(ma))))
  if (any(duplicated(all_scalar_names))) {
    dupes <- unique(all_scalar_names[duplicated(all_scalar_names)])
    stop(
      "Scalar name collision across ModelArrays: ",
      paste(dupes, collapse = ", "),
      ". Each scalar must have a unique name."
    )
  }

  # Inner join phenotypes on merge_on columns
  merged_phen <- phenotypes_list[[1]]
  for (i in 2:length(phenotypes_list)) {
    merged_phen <- merge(merged_phen, phenotypes_list[[i]],
      by = merge_on, suffixes = c("", paste0(".", i))
    )
  }
  if (nrow(merged_phen) == 0) {
    stop("Inner join on merge_on columns produced zero rows. No shared subjects found.")
  }

  # Create a unified source identifier from the merge_on columns.
  # All scalars in the combined ModelArray will share this identifier,
  # which allows cross-scalar formulas in analyseOneElement.* to work.
  if (length(merge_on) == 1) {
    unified_source_id <- as.character(merged_phen[[merge_on]])
  } else {
    unified_source_id <- apply(
      merged_phen[, merge_on, drop = FALSE], 1,
      function(row) paste(row, collapse = "_")
    )
  }
  if (any(duplicated(unified_source_id))) {
    stop(
      "merge_on columns do not uniquely identify subjects. ",
      "Ensure the combination of merge_on columns is unique per subject."
    )
  }

  # Build combined scalars, sources, and path
  combined_scalars <- list()
  combined_sources <- list()
  combined_paths <- character()

  for (i in seq_along(modelarrays)) {
    ma <- modelarrays[[i]]
    ma_scalar_names <- names(scalars(ma))

    # Determine which source_file column corresponds to this ModelArray
    if (i == 1) {
      sf_col <- "source_file"
    } else {
      sf_col <- paste0("source_file.", i)
    }

    # Get the source files from the merged phenotypes for this ModelArray
    original_sources <- merged_phen[[sf_col]]

    for (sn in ma_scalar_names) {
      # Map original sources to h5 column indices
      ma_sources <- sources(ma)[[sn]]
      col_idx <- match(original_sources, ma_sources)
      if (any(is.na(col_idx))) {
        missing <- original_sources[is.na(col_idx)]
        stop(
          "Source files from merged phenotypes not found in ModelArray for scalar '",
          sn, "': ", paste(head(missing, 3), collapse = ", ")
        )
      }

      # Subset the DelayedArray columns to match the merged subject order
      # Use unified_source_id as the column names for all scalars
      combined_scalars[[sn]] <- scalars(ma)[[sn]][, col_idx, drop = FALSE]
      colnames(combined_scalars[[sn]]) <- unified_source_id
      combined_sources[[sn]] <- unified_source_id
      combined_paths[sn] <- ma@path[1]
    }
  }

  # Check element correspondence across all scalars
  n_elements <- sapply(combined_scalars, nrow)
  if (length(unique(n_elements)) > 1) {
    detail <- paste(names(n_elements), n_elements, sep = "=", collapse = ", ")
    stop("Element count mismatch across ModelArrays: ", detail)
  }

  # Check element metadata if available
  metadata_list <- lapply(modelarrays, function(ma) {
    tryCatch(elementMetadata(ma), error = function(e) NULL)
  })
  non_null_meta <- Filter(Negate(is.null), metadata_list)
  if (length(non_null_meta) >= 2) {
    ref_meta <- non_null_meta[[1]]
    for (j in 2:length(non_null_meta)) {
      meta_differs <- !identical(dim(ref_meta), dim(non_null_meta[[j]])) ||
        !all(ref_meta == non_null_meta[[j]], na.rm = TRUE)
      if (meta_differs) {
        warning(
          "Element metadata differs between ModelArrays. ",
          "Ensure elements correspond to the same spatial locations.",
          call. = FALSE
        )
        break
      }
    }
  } else if (length(non_null_meta) < length(modelarrays) && length(non_null_meta) > 0) {
    warning(
      "Element metadata not available in all ModelArrays. ",
      "Cannot verify element correspondence across files.",
      call. = FALSE
    )
  }

  # Rename original source_file columns to be scalar-specific,
  # and add unified source_file column for use with analysis functions
  first_scalar_per_ma <- sapply(modelarrays, function(ma) names(scalars(ma))[1])
  for (i in seq_along(modelarrays)) {
    if (i == 1) {
      old_name <- "source_file"
    } else {
      old_name <- paste0("source_file.", i)
    }
    new_name <- paste0("source_file.", first_scalar_per_ma[i])
    idx <- which(colnames(merged_phen) == old_name)
    if (length(idx) == 1) {
      colnames(merged_phen)[idx] <- new_name
    }
  }
  # The unified source_file column matches all scalars' sources
  merged_phen[["source_file"]] <- unified_source_id

  # Construct the combined ModelArray
  combined_ma <- new("ModelArray",
    scalars = combined_scalars,
    sources = combined_sources,
    results = list(),
    path = combined_paths
  )

  list(data = combined_ma, phenotypes = merged_phen)
}
