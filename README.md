# Kinetochore spot intensity

Fiji/ImageJ macros that measure the intensity of a protein of interest at individual kinetochores, normalized to a kinetochore reference marker (e.g., CENP-A).

Used in Skendo et al. (Cell Reports) to quantify Nnf1R and Hec1/NDC80 at kinetochores relative to CENP-A in hTERT-RPE1 Centrin1-GFP/GFP-CENP-A cells (Figures S3C and S3D).

| Step | Macro | What it does |
|------|-------|--------------|
| 1 | `Kinetochore_Spot_Intensity_Extractor.ijm` | Detects kinetochores on the reference channel, measures both channels at each kinetochore, subtracts background and computes the protein/reference ratio. One CSV per cell. |
| 2 | `KT_Spot_PerCell_Summary.ijm` | Averages all kinetochores within each cell into one row per cell, for two conditions. |

## Step 1: per-kinetochore measurement

1. Open a multi-channel sum projection in Fiji.
2. Run `Kinetochore_Spot_Intensity_Extractor.ijm` (**Plugins > Macros > Run**).
3. Assign the protein-of-interest channel and the poles + kinetochores reference channel. Set the Find Maxima prominence (default 500) and the spot radius (default 3 pixels). In Skendo et al., the radius was 3 pixels and the prominence was 500, lowered to 250 for images where too few kinetochores were detected.
4. Find Maxima runs on the reference channel, so kinetochore detection does not depend on the protein of interest. The candidate points are added to the ROI Manager: delete the spindle pole points, then click **OK**.
5. Draw a rectangular background ROI outside the spindle and click **OK**.
6. Choose an output folder. `<imagename>_KTspots.csv` is saved with one row per kinetochore:

`Image, SpotIndex, X_px, Y_px, Protein_raw, Protein_bgsub, Ref_raw, Ref_bgsub, Ratio_bgsub`

`Ratio_bgsub` = `Protein_bgsub` / `Ref_bgsub`. It is left blank when `Ref_bgsub` is close to zero.

## Step 2: per-cell summary

1. Put all `*_KTspots.csv` files (all conditions and replicates) in one folder.
2. Run `KT_Spot_PerCell_Summary.ijm` and select the folder.
3. Enter a name and one or more filename keywords for each of the two conditions (e.g., `siCTRL` / `ctrl` and `siMIS12` / `mis12`).
4. `percell_KTspot_summary.csv` (and `.xls`) is saved in the same folder:

`Filename, Condition, N_kinetochores, MeanProtein_bgsub, MeanRef_bgsub, MeanRatio_bgsub`

Kinetochores within one cell are not independent, so the cell is the unit for statistics. Per-cell values were then averaged per experiment and normalized to siCTRL of the same experiment.

## Requirements

[Fiji/ImageJ](https://fiji.sc)

## Author

Kristjana Skendo, Meraldi Lab, University of Geneva
