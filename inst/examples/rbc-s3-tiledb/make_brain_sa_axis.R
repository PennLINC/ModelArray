#!/usr/bin/env Rscript
# Build the S-A axis brain images that ship with the example: lateral views of the left
# inflated fsLR-164k surface coloured by the Glasser sensorimotor-association (S-A) axis,
# in the colour scheme of Sydnor et al. 2021 (Neuron), Figure 2A -- gold (sensorimotor)
# through white to purple (association), diverging at rank 90.
#
#   make_brain_sa_axis.R <L inflated 164k surf.gii> <glasser360_fsLR164k_labels.rds> <out dir> [wb_command]
#
# Writes three PNGs into <out dir>:
#   brain_sa_axis.png  areas 3b and 46 outlined and labelled (the main locator)
#   brain_sa_3b.png    area 3b outlined only (for the 3b panel of the vertex figure)
#   brain_sa_46.png    area 46 outlined only (for the 46 panel of the vertex figure)
#
# The 164k inflated surface is the one the example builds from TemplateFlow's 32k surfaces
# (see the report's cifti display step, or make_glasser_fsLR164k_labels.sh). This is
# provenance for files that ship with the example; users do not run it. Needs wb_command
# (to convert the GIFTI to ASCII), xml2, ggplot2 and scales.
#
# S-A axis ranks: Sydnor VJ et al. 2021, Neuron (PennLINC/S-A_ArchetypalAxis). Cite Sydnor
# et al. 2021 if you use them. The surface is rendered as shaded, back-face-culled filled
# triangles; parcel outlines are the mesh edges between in- and out-of-parcel triangles.
suppressMessages({library(xml2); library(ggplot2); library(scales)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: make_brain_sa_axis.R L_inflated.surf.gii labels.rds out_dir [wb_command]")
surf_gii <- args[1]; labels_rds <- args[2]; out_dir <- args[3]
wb <- if (length(args) >= 4) args[4] else "wb_command"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ascii <- tempfile(fileext = ".surf.gii")
stopifnot(system2(wb, c("-gifti-convert", "ASCII", surf_gii, ascii)) == 0)
doc <- read_xml(ascii)
V <- matrix(as.numeric(scan(text = xml_text(xml_find_first(doc, ".//DataArray[@Intent='NIFTI_INTENT_POINTSET']//Data")), quiet = TRUE)), ncol = 3, byrow = TRUE)
Tt <- matrix(as.integer(scan(text = xml_text(xml_find_first(doc, ".//DataArray[@Intent='NIFTI_INTENT_TRIANGLE']//Data")), quiet = TRUE)), ncol = 3, byrow = TRUE) + 1L
labels <- readRDS(labels_rds)
keyL <- labels$key[seq_len(nrow(V))]
stopifnot(nrow(V) == 163842L)

SA_COLS <- c("goldenrod1", "white", "#6f1282"); SA_MID <- 90; SA_RANGE <- c(1, 180)   # Sydnor 2021 Fig 2A
rank_col <- function(r) gradient_n_pal(SA_COLS)(rescale_mid(r, to = c(0, 1), from = SA_RANGE, mid = SA_MID))
WALL <- "grey92"
rankv <- labels$table$sa_rank[match(keyL, labels$table$key)]
sacol <- ifelse(is.na(rankv), WALL, rank_col(rankv))

blend <- function(hex, f) { m <- grDevices::col2rgb(hex) / 255; m <- sweep(m, 2, f, `*`); m[m > 1] <- 1; m[m < 0] <- 0; grDevices::rgb(m[1,], m[2,], m[3,]) }
# left lateral view: anterior (+Y) to the left, superior (+Z) up; keep lateral-facing triangles
xi <- -V[, 2]; yi <- V[, 3]
a <- V[Tt[,1],]; b <- V[Tt[,2],]; c <- V[Tt[,3],]
n <- cbind((b[,2]-a[,2])*(c[,3]-a[,3])-(b[,3]-a[,3])*(c[,2]-a[,2]),
           (b[,3]-a[,3])*(c[,1]-a[,1])-(b[,1]-a[,1])*(c[,3]-a[,3]),
           (b[,1]-a[,1])*(c[,2]-a[,2])-(b[,2]-a[,2])*(c[,1]-a[,1]))
nn <- sqrt(rowSums(n^2)); nn[nn == 0] <- 1; n <- n / nn
ff <- n[, 1] < 0; Tf <- Tt[ff, , drop = FALSE]; nf <- n[ff, , drop = FALSE]
light <- c(-0.5, 0.25, 0.85); light <- light / sqrt(sum(light^2))
shade <- pmax(0.6, 0.78 + 0.22 * as.vector(nf %*% light))
tcol <- blend(sacol[Tf[, 1]], shade)
depth <- (V[Tf[,1],1] + V[Tf[,2],1] + V[Tf[,3],1]) / 3
ord <- order(depth, decreasing = TRUE); Tf <- Tf[ord, ]; tcol <- tcol[ord]; m <- nrow(Tf)
poly <- data.frame(x = c(xi[Tf[,1]], xi[Tf[,2]], xi[Tf[,3]]), y = c(yi[Tf[,1]], yi[Tf[,2]], yi[Tf[,3]]),
                   grp = rep(seq_len(m), 3), fill = rep(tcol, 3))
poly <- poly[order(poly$grp), ]
base_plot <- ggplot(poly, aes(x, y, group = grp, fill = fill)) + geom_polygon(colour = NA) +
  scale_fill_identity() + coord_fixed() + theme_void()

outline_df <- function(nm) {
  inpar <- keyL == labels$table$key[labels$table$name == nm]
  tri_in <- rowSums(matrix(inpar[Tf], ncol = 3)) >= 2L
  e <- rbind(Tf[,c(1,2)], Tf[,c(2,3)], Tf[,c(1,3)]); e <- cbind(pmin(e[,1],e[,2]), pmax(e[,1],e[,2]))
  k <- paste(e[,1], e[,2])
  keep <- ave(rep(1L, nrow(e)), k, FUN = sum) == 2L & ave(as.integer(rep(tri_in, 3)), k, FUN = sum) == 1L & !duplicated(k)
  e <- e[keep, , drop = FALSE]
  data.frame(x = xi[e[,1]], y = yi[e[,1]], xend = xi[e[,2]], yend = yi[e[,2]])
}
centroid <- function(nm) { v <- which(keyL == labels$table$key[labels$table$name == nm]); data.frame(x = mean(xi[v]), y = mean(yi[v])) }

# region name -> (glasser parcel, label text, label nudge)
PARCELS <- list("3b" = list(roi = "L_3b_ROI", dx = -14, dy = 16),
                "46" = list(roi = "L_46_ROI", dx = -2,  dy = 18))
render <- function(regions, label) {
  g <- base_plot
  for (rg in regions) g <- g + geom_segment(data = outline_df(PARCELS[[rg]]$roi),
    aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "black", linewidth = 0.45)
  if (label) {
    labs <- do.call(rbind, lapply(regions, function(rg) cbind(centroid(PARCELS[[rg]]$roi), lab = rg,
      dx = PARCELS[[rg]]$dx, dy = PARCELS[[rg]]$dy)))
    g <- g + geom_text(data = labs, aes(x + dx, y + dy, label = lab), inherit.aes = FALSE, colour = "black", fontface = 2, size = 5.4)
  }
  g
}
save_png <- function(g, file) {
  ggsave(file.path(out_dir, file), g, width = 3.4, height = 2.8, dpi = 200, bg = "transparent")
  message(sprintf("wrote %s (%.0f KB)", file, file.size(file.path(out_dir, file)) / 1024))
}
save_png(render(c("3b", "46"), label = TRUE),  "brain_sa_axis.png")
save_png(render("3b", label = FALSE),          "brain_sa_3b.png")
save_png(render("46", label = FALSE),          "brain_sa_46.png")
