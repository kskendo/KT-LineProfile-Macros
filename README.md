# KT-LineProfile-Macros

**Repository:** https://github.com/kskendo/KT-LineProfile-Macros

Fiji/ImageJ macros for kinetochore-to-kinetochore (KT–KT) fluorescence intensity line profile analysis across bridging fiber overlaps.

Developed for:
> Skendo et al. (2025) *Title to be confirmed upon publication.* Cell Reports.

---

## Overview

These two macros work together as a pipeline:

| Step | Macro | What it does |
|------|-------|--------------|
| 1 | `Line-profiles_2KTpair_final.ijm` | Extracts intensity profiles along the KT–KT axis for 2 kinetochore pairs per cell, across 3 channels. Outputs per-pair and per-cell average profiles as a CSV. |
| 2 | `averaging-line-profiles.ijm` | Takes all per-image CSV files from Step 1 and computes the grand mean ± SEM across all cells. |

---

## Requirements

- [Fiji/ImageJ](https://fiji.sc) (any recent version)
- Multi-channel z-stack images (`.tif`, `.tiff`, or `.dv` format)
- Images should be deconvolved and contain at least 3 channels:
  - Channel 1: kinetochore marker (e.g., CENP-A/GFP)
  - Channel 2: protein of interest 1 (e.g., HURP or KID)
  - Channel 3: protein of interest 2 (e.g., PRC1 or α-tubulin)

---

## Usage

### Step 1 — Extract per-cell profiles (`Line-profiles_2KTpair_final.ijm`)

1. Open a multi-channel z-stack in Fiji.
2. Run the macro via **Plugins > Macros > Run**.
3. Enter the following when prompted:
   - Channel numbers for kinetochore marker and two proteins
   - Number of bins for profile resampling (default: 101)
   - Number of cells to analyze in this image
   - Line width in pixels
4. For each cell, for each of the 2 KT pairs:
   - Navigate to the z-slice where both kinetochores are clearly visible
   - Draw a straight line from one kinetochore to the other using the line tool
   - Click **OK**
5. Choose an output folder. A CSV named `<imagename>_profiles.csv` is saved automatically.

**Output CSV structure** (`<imagename>_profiles.csv`):

| Column | Description |
|--------|-------------|
| `image` | Image filename (without extension) |
| `cell` | Cell number |
| `pair_type` | `pair1`, `pair2`, or `average` (mean of the two pairs for that cell) |
| `position_0to1` | Normalized position along KT–KT axis (0 = first KT, 1 = second KT) |
| `centered_pos` | Same, centered at 0 (range: −0.5 to 0.5) |
| `ch1_raw` / `ch2_raw` / `ch3_raw` | Raw intensity at each bin |
| `ch1_norm` / `ch2_norm` / `ch3_norm` | Min–max normalized intensity (0–1) per profile |

**Example rows:**
```
image,cell,pair_type,position_0to1,centered_pos,ch1_raw,ch2_raw,ch3_raw,ch1_norm,ch2_norm,ch3_norm
KID_PRC1_1_01_03_R3D_D3D,1,pair1,0.000000,-0.500000,1277.312500,362.250000,470.062500,0.489166,0.000000,0.000000
KID_PRC1_1_01_03_R3D_D3D,1,pair1,0.010000,-0.490000,1361.484230,366.250830,472.921981,0.543460,0.012447,0.013679
KID_PRC1_1_01_03_R3D_D3D,1,average,0.000000,-0.500000,1160.406250,424.875000,499.093750,0.506111,0.000000,0.000000
```

Each image with 1 cell and 101 bins produces 303 rows (101 × 3 pair types).

---

### Step 2 — Compute grand mean ± SEM (`averaging-line-profiles.ijm`)

1. Place all CSV files from Step 1 for **one condition** into a single folder.
2. Run the macro via **Plugins > Macros > Run**.
3. Select the folder containing the CSV files.
4. Output is saved in the same folder as:
   `grand_average_per_cell_normalized_profiles_with_SEM.csv`

Only rows with `pair_type = average` are used — one data point per cell per bin — from which the grand mean and SEM are calculated across all cells.

**Output CSV structure:**

| Column | Description |
|--------|-------------|
| `position_0to1` | Normalized position (0–1) |
| `centered_pos` | Centered position (−0.5 to 0.5) |
| `ch1_mean` / `ch2_mean` / `ch3_mean` | Grand mean normalized intensity across all cells |
| `ch1_sem` / `ch2_sem` / `ch3_sem` | SEM across all cells (Bessel's correction, n−1) |
| `n_cells` | Number of cells used |

**Example rows:**
```
position_0to1,centered_pos,ch1_mean,ch2_mean,ch3_mean,ch1_sem,ch2_sem,ch3_sem,n_cells
0.000000,-0.500000,0.481037,0.407219,0.195786,0.050866,0.071599,0.047211,30
0.010000,-0.490000,0.515510,0.411445,0.211127,0.048688,0.069247,0.046851,30
0.500000,0.000000,0.025917,0.561754,0.621923,0.005701,0.054011,0.047446,30
```

---

## Workflow

```
Image 1 (.dv/.tif) ──► Step 1 ──► image1_profiles.csv ─┐
Image 2 (.dv/.tif) ──► Step 1 ──► image2_profiles.csv ─┼──► Step 2 ──► grand_average_with_SEM.csv
Image 3 (.dv/.tif) ──► Step 1 ──► image3_profiles.csv ─┘
         (one folder per condition)
```

Run Step 2 separately for each experimental condition (e.g., siControl and siHURP in separate folders).

---

## Notes

- The number of bins (default: 101) must be consistent across all images fed into Step 2.
- Line width (in pixels) sets the averaging width for the intensity profile — adjust based on your image resolution and structure width.
- Profiles are aligned to the KT–KT axis only; no pole-to-pole alignment is applied.
- The `centered_pos` column is recommended for plotting, as it centers the profile at the midpoint between the two kinetochores (position 0).
- For downstream plotting, use the `*_norm` columns from Step 1 and `*_mean`/`*_sem` columns from Step 2.

---

## Example data

Two example output files are provided in the `example_data/` folder:
- `KID_PRC1_1_01_03_R3D_D3D_profiles.csv` — example output from Step 1 (1 cell, 3 channels: CENP-A, KID, PRC1)
- `grand_average_per_cell_normalized_profiles_with_SEM.csv` — example output from Step 2 (30 cells)

---

## Citation

If you use these macros, please cite:

> Skendo K et al. (2025) Title to be confirmed. *Cell Reports*. DOI: to be added upon publication.

---

## License

MIT License — free to use and modify with attribution. See `LICENSE` for details.

---

## Contact

Kristjana Skendo  
Meraldi Lab, Department of Cell Physiology and Metabolism  
University of Geneva  
patrick.meraldi@unige.ch
