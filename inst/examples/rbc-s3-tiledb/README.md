# RBC cortical thickness from S3: a runnable ModelArray example

`rbc_s3_tiledb_example.Rmd` reads the RBC FreeSurferPost TileDB store straight from
its public S3 bucket, fits vertex-wise GAMs of cortical thickness on age in a handful
of Glasser parcels, and writes the results as CIFTI files you can open in Connectome
Workbench. `vignette("rbc-s3-tiledb")` explains how to run it and what to expect.

## Files

| File | What it is |
|---|---|
| `rbc_s3_tiledb_example.Rmd` | The report. Render it with `rmarkdown::render()`; see the vignette for the exact call and parameters. |
| `glasser360_fsLR164k_labels.rds` | One Glasser (HCP-MMP1.0) parcel key per column of the TileDB store (327,684 fsLR-164k vertices, left hemisphere then right, medial wall = 0), with the label table and each parcel's sensorimotor-association (S-A) axis rank. 26 KB. |
| `sa_axis_glasser360_ranks.csv` | S-A axis average ranks for the 180 Glasser regions (`final.rank`: 1 = sensorimotor pole, 180 = association pole). |
| `make_glasser_fsLR164k_labels.sh`, `make_glasser_fsLR164k_labels.R` | How the `.rds` was made: the 32k Glasser dlabel resampled to 164k with Connectome Workbench and TemplateFlow's fsLR spheres. Kept for provenance; you don't need to run them. |
| `brain_sa_axis.png`, `brain_sa_3b.png`, `brain_sa_46.png` | Lateral views of the left inflated cortex coloured by the S-A axis in the colour scheme of Sydnor et al. 2021 Figure 2A (gold = sensorimotor, purple = association). `brain_sa_axis` outlines and labels both areas 3b and 46; the other two outline one parcel each. The figures use them as locators. |
| `make_brain_sa_axis.R` | How the three brain PNGs were made (filled-mesh render of the 164k inflated surface). Provenance; you don't need to run it. |

## Sources and attribution

- **Parcellation:** Glasser MF et al. (2016). A multi-modal parcellation of human
  cerebral cortex. *Nature* 536, 171–178. The 32k dlabel is the copy distributed with
  [XCP-D](https://github.com/PennLINC/xcp_d) (`atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii`);
  the HCP-MMP1.0 parcellation is released under the HCP open-access data use terms.
- **S-A axis ranks:** Sydnor VJ et al. (2021). Neurodevelopment of the association
  cortices: patterns, mechanisms, and implications for psychopathology. *Neuron* 109,
  2820–2846. File from [PennLINC/S-A_ArchetypalAxis](https://github.com/PennLINC/S-A_ArchetypalAxis)
  (`Glasser360_MMP/Sensorimotor_Association_Axis_AverageRanks.csv`). Please cite
  Sydnor et al. 2021 if you use the ranks. The S-A axis brain (`brain_sa_axis.png`)
  reuses that paper's Figure 2A colour scheme (gold → white → purple, diverging at
  rank 90).
- **Surfaces for resampling and display:** [TemplateFlow](https://www.templateflow.org/)
  `tpl-fsLR` (spheres, vertex-area maps and 32k display surfaces), downloaded by the
  scripts when needed.
- **Data:** Reproducible Brain Charts (RBC) FreeSurferPost release, published as a
  TileDB store at `s3://fcp-indi/data/Projects/RBC/RBC_FreeSurfer-Post_ModelArray/`.
  See the RBC documentation at <https://reprobrainchart.github.io/> for the data
  release and its terms.
