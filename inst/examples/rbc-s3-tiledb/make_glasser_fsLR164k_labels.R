#!/usr/bin/env Rscript
# Second half of make_glasser_fsLR164k_labels.sh: turn the text export of the 164k
# Glasser dlabel into glasser360_fsLR164k_labels.rds, with checks.
#
#   Rscript make_glasser_fsLR164k_labels.R <keys.txt> <labeltable.txt> <sa_ranks.csv> <out.rds> <input dlabel> <wb version>
#
# The .rds holds a list:
#   key        integer vector, one entry per TileDB column (327,684): 0 = medial wall / unlabelled
#   table      data.frame key, name, hemi, region, sa_rank: one row per parcel (360)
#   provenance list: inputs, method, date, software
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: make_glasser_fsLR164k_labels.R keys.txt labeltable.txt sa_ranks.csv out.rds input.dlabel.nii [wb_version]")
keys_txt <- args[1]; table_txt <- args[2]; sa_csv <- args[3]; out_rds <- args[4]; input_dlabel <- args[5]
wb_version <- if (length(args) >= 6) args[6] else "unknown"

key <- as.integer(scan(keys_txt, what = numeric(), quiet = TRUE))
n_hemi <- 163842L
stopifnot(length(key) == 2L * n_hemi)

# -cifti-label-export-table writes pairs of lines: the name, then "key R G B A"
lines <- readLines(table_txt)
lines <- lines[nzchar(trimws(lines))]
stopifnot(length(lines) %% 2 == 0)
names_ <- lines[seq(1, length(lines), by = 2)]
keys_ <- as.integer(vapply(strsplit(trimws(lines[seq(2, length(lines), by = 2)]), "[[:space:]]+"), `[`, "", 1))
table <- data.frame(key = keys_, name = names_, stringsAsFactors = FALSE)
table <- table[table$key != 0L, ]                       # drop the "???" background label
table$hemi <- substr(table$name, 1, 1)
table$region <- sub("^[LR]_(.*)_ROI$", "\\1", table$name)
stopifnot(nrow(table) == 360L, all(table$hemi %in% c("L", "R")), length(unique(table$region)) == 180L)

# Sensorimotor-association axis rank per region (Sydnor et al. 2021), same for both hemispheres
sa <- read.csv(sa_csv, stringsAsFactors = FALSE)
table$sa_rank <- sa$final.rank[match(table$region, sa$region)]
missing_rank <- table$region[is.na(table$sa_rank)]
if (length(missing_rank) > 0) stop("No S-A rank for: ", paste(unique(missing_rank), collapse = ", "))

# Alignment with the TileDB column order: left vertices first, then right
left_keys <- table$key[table$hemi == "L"]
stopifnot(
  all(key[seq_len(n_hemi)] %in% c(0L, left_keys)),
  all(key[(n_hemi + 1):(2 * n_hemi)] %in% c(0L, table$key[table$hemi == "R"])),
  all(table$key %in% key)
)
table <- table[order(table$key), ]
rownames(table) <- NULL

counts <- table(factor(key, levels = c(0L, table$key)))
message(sprintf("medial wall / unlabelled: %d of %d vertices; parcel sizes %d-%d vertices (median %d)",
                counts[["0"]], length(key), min(counts[-1]), max(counts[-1]), as.integer(median(counts[-1]))))

out <- list(
  key = key,
  table = table,
  provenance = list(
    description = "Glasser et al. 2016 HCP-MMP1.0 parcellation on fsLR 164k, one key per grayordinate of the RBC FreeSurferPost dscalars (CORTEX_LEFT 0-163841, then CORTEX_RIGHT; medial wall = 0)",
    input_dlabel = basename(input_dlabel),
    input_md5 = unname(tools::md5sum(input_dlabel)),
    method = "wb_command -cifti-separate; -label-resample ADAP_BARY_AREA with TemplateFlow tpl-fsLR 32k/164k spheres and desc-vaavg_midthickness vertex areas; -cifti-create-label; -cifti-convert -to-text",
    sa_axis = "Sensorimotor_Association_Axis_AverageRanks.csv, PennLINC/S-A_ArchetypalAxis (Sydnor et al. 2021, Neuron), column final.rank",
    workbench = wb_version,
    date = format(Sys.Date())
  )
)
saveRDS(out, out_rds, compress = "xz")
message(sprintf("wrote %s (%.0f KB)", out_rds, file.size(out_rds) / 1024))
