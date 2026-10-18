#!/usr/bin/env bash
# Build glasser360_fsLR164k_labels.rds: one Glasser (HCP-MMP1.0) parcel key for each of
# the 327,684 columns of the RBC FreeSurferPost TileDB store (fsLR 164k, all vertices,
# left hemisphere then right, medial wall = 0).
#
# This is provenance for a file that ships with the example; users of the example do
# not run it. It needs Connectome Workbench (wb_command), curl and Rscript.
#
#   make_glasser_fsLR164k_labels.sh <atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii> <workdir>
#
# The 32k Glasser dlabel is the one distributed with XCP-D
# (xcp_d/data/atlases/atlas-Glasser/). Spheres and vertex-area maps come from
# TemplateFlow's tpl-fsLR; the grayordinate layout is checked against one release
# dscalar, so the text export is in the same order as the TileDB columns.
set -euo pipefail
IN=${1:?usage: make_glasser_fsLR164k_labels.sh <glasser 32k dlabel> <workdir>}
W=${2:?usage: make_glasser_fsLR164k_labels.sh <glasser 32k dlabel> <workdir>}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WB=${WB_COMMAND:-wb_command}
TF=https://templateflow.s3.amazonaws.com/tpl-fsLR
TEMPLATE_URL=https://fcp-indi.s3.amazonaws.com/data/Projects/RBC/CCNP_FreeSurfer-Post/SHA256E-s3385504--99c174365d91c22508a0de725b40e19f8e6b572c6efd88dfcd681e392e187267.nii
TEMPLATE_SHA=99c174365d91c22508a0de725b40e19f8e6b572c6efd88dfcd681e392e187267

mkdir -p "$W"
for h in L R; do
  for d in 32k 164k; do
    for f in "tpl-fsLR_hemi-${h}_den-${d}_sphere.surf.gii" "tpl-fsLR_hemi-${h}_den-${d}_desc-vaavg_midthickness.shape.gii"; do
      [ -s "$W/$f" ] || curl -sSfL -o "$W/$f" "$TF/$f"
    done
  done
done
[ -s "$W/template_fsLR164k.dscalar.nii" ] || curl -sSfL -o "$W/template_fsLR164k.dscalar.nii" "$TEMPLATE_URL"
echo "$TEMPLATE_SHA  $W/template_fsLR164k.dscalar.nii" | shasum -a 256 -c - >/dev/null

# 32k dlabel -> one label.gii per hemisphere -> 164k with area-corrected barycentric resampling
"$WB" -cifti-separate "$IN" COLUMN -label CORTEX_LEFT "$W/L_32k.label.gii" -label CORTEX_RIGHT "$W/R_32k.label.gii"
for h in L R; do
  "$WB" -label-resample "$W/${h}_32k.label.gii" \
    "$W/tpl-fsLR_hemi-${h}_den-32k_sphere.surf.gii" "$W/tpl-fsLR_hemi-${h}_den-164k_sphere.surf.gii" \
    ADAP_BARY_AREA "$W/${h}_164k.label.gii" \
    -area-metrics "$W/tpl-fsLR_hemi-${h}_den-32k_desc-vaavg_midthickness.shape.gii" \
                  "$W/tpl-fsLR_hemi-${h}_den-164k_desc-vaavg_midthickness.shape.gii"
done
"$WB" -cifti-create-label "$W/glasser_fsLR164k.dlabel.nii" \
  -left-label "$W/L_164k.label.gii" -right-label "$W/R_164k.label.gii"

# Both files must list all 163,842 vertices per hemisphere, left before right
for f in glasser_fsLR164k.dlabel.nii template_fsLR164k.dscalar.nii; do
  echo "$f: $("$WB" -file-information "$W/$f" | grep -E 'Number of Rows|out of' | tr -s ' ' | paste -sd' ' -)"
done

"$WB" -cifti-convert -to-text "$W/glasser_fsLR164k.dlabel.nii" "$W/glasser_fsLR164k_keys.txt"
"$WB" -cifti-label-export-table "$W/glasser_fsLR164k.dlabel.nii" 1 "$W/glasser_fsLR164k_labeltable.txt"
Rscript "$HERE/make_glasser_fsLR164k_labels.R" \
  "$W/glasser_fsLR164k_keys.txt" "$W/glasser_fsLR164k_labeltable.txt" \
  "$HERE/sa_axis_glasser360_ranks.csv" "$HERE/glasser360_fsLR164k_labels.rds" \
  "$IN" "$("$WB" -version | head -1)"
