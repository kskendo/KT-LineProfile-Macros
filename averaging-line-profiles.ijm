// ============================================================
// Macro: Grand Average Line Profile Calculator (with SEM)
// ============================================================
// Description:
//   Takes a folder of CSV files produced by the companion macro
//   "Line-profiles_2KTpair_final.ijm" and computes the grand
//   mean and SEM of the per-cell averaged normalized profiles
//   across all cells and all CSV files in the folder.
//
//   Only rows with pair_type = "average" are used (i.e., the
//   per-cell average of the two KT pairs). This gives one data
//   point per cell per bin, from which the grand mean and SEM
//   across all cells are calculated.
//
// Usage:
//   1. Run "Line-profiles_2KTpair_final.ijm" on all your images
//      to generate per-image CSV files. Place them all in one folder.
//   2. Run this macro via Plugins > Macros > Run.
//   3. Select the folder containing the CSV files.
//   4. Output is saved in the same folder as:
//      "grand_average_per_cell_normalized_profiles_with_SEM.csv"
//
// Input:
//   - Folder of CSV files from "Line-profiles_2KTpair_final.ijm"
//   - Each CSV must contain columns:
//     image, cell, pair_type, position_0to1, centered_pos,
//     ch1_raw, ch2_raw, ch3_raw, ch1_norm, ch2_norm, ch3_norm
//   - Only rows where pair_type == "average" are processed
//
// Output:
//   - grand_average_per_cell_normalized_profiles_with_SEM.csv with columns:
//     position_0to1, centered_pos,
//     ch1_mean, ch2_mean, ch3_mean,
//     ch1_sem, ch2_sem, ch3_sem,
//     n_cells
//
// Notes:
//   - SEM is calculated as: SD / sqrt(n_cells)
//   - All CSV files in the selected folder are processed together,
//     so one folder = one experimental condition
//   - Bin count is detected automatically from the first cell encountered
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// Publication: Skendo et al., 2025 (Cell Reports, title to be confirmed)
// Repository: https://github.com/kskendo/KT-LineProfile-Macros
// ============================================================

macro "Average normalized average profiles per cell + SEM" {

    inputDir = getDirectory("Choose folder with CSV files");
    if (inputDir == "") exit("No input folder selected.");

    fileList = getFileList(inputDir);
    if (fileList.length == 0) exit("No files found in the folder.");

    outFile = inputDir + "grand_average_per_cell_normalized_profiles_with_SEM.csv";

    // --- Accumulators ---
    initialized = 0;
    nCells = 0;

    // --- Loop over all CSV files in the folder ---
    for (f = 0; f < fileList.length; f++) {
        name = fileList[f];
        if (!endsWith(name, ".csv")) continue;

        path = inputDir + name;
        text = File.openAsString(path);
        lines = split(text, "\n");

        currentCellKey = "";
        binIndex = 0;

        // --- Loop over rows in this CSV ---
        for (i = 1; i < lines.length; i++) {
            line = trim(lines[i]);
            if (line == "") continue;

            cols = split(line, ",");
            if (cols.length < 11) continue;

            imageName = cols[0];
            cellID    = cols[1];
            pairType  = cols[2];

            // Only use per-cell average rows
            if (pairType != "average") continue;

            // Unique key per cell (image + cell ID)
            cellKey = imageName + "__cell" + cellID;

            // Detect start of a new cell
            if (cellKey != currentCellKey) {
                currentCellKey = cellKey;
                binIndex = 0;
                nCells++;
            }

            // Parse values (use normalized columns: ch1_norm=col8, ch2_norm=col9, ch3_norm=col10)
            posVal      = parseFloat(cols[3]);
            centeredVal = parseFloat(cols[4]);
            ch1 = parseFloat(cols[8]);
            ch2 = parseFloat(cols[9]);
            ch3 = parseFloat(cols[10]);

            // Initialize accumulator arrays on first row
            if (initialized == 0) {
                tempPos      = newArray(0);
                tempCentered = newArray(0);
                tempSum1     = newArray(0); tempSum2 = newArray(0); tempSum3 = newArray(0);
                tempSumsq1   = newArray(0); tempSumsq2 = newArray(0); tempSumsq3 = newArray(0);
                initialized  = 1;
            }

            // Grow accumulator arrays when a new bin is encountered
            if (binIndex >= lengthOf(tempPos)) {
                tempPos      = appendValue(tempPos, posVal);
                tempCentered = appendValue(tempCentered, centeredVal);
                tempSum1     = appendValue(tempSum1, 0);
                tempSum2     = appendValue(tempSum2, 0);
                tempSum3     = appendValue(tempSum3, 0);
                tempSumsq1   = appendValue(tempSumsq1, 0);
                tempSumsq2   = appendValue(tempSumsq2, 0);
                tempSumsq3   = appendValue(tempSumsq3, 0);
            }

            // Warn if bin positions are inconsistent across files
            if (abs(tempPos[binIndex] - posVal) > 0.000001) {
                print("Warning: position mismatch in file " + name + ", cell " + cellID + ", bin " + binIndex);
            }

            // Accumulate sum and sum-of-squares for mean/SEM calculation
            tempSum1[binIndex] += ch1;   tempSum2[binIndex] += ch2;   tempSum3[binIndex] += ch3;
            tempSumsq1[binIndex] += ch1 * ch1;
            tempSumsq2[binIndex] += ch2 * ch2;
            tempSumsq3[binIndex] += ch3 * ch3;

            binIndex++;
        }
    }

    if (initialized == 0 || nCells == 0)
        exit("No rows with pair_type = 'average' were found. Check your CSV files.");

    nBins = lengthOf(tempPos);

    // --- Build output CSV ---
    output = "position_0to1,centered_pos,ch1_mean,ch2_mean,ch3_mean,ch1_sem,ch2_sem,ch3_sem,n_cells\n";

    for (j = 0; j < nBins; j++) {
        // Grand mean per bin
        mean1 = tempSum1[j] / nCells;
        mean2 = tempSum2[j] / nCells;
        mean3 = tempSum3[j] / nCells;

        // SEM = sqrt(variance / n), using Bessel's correction (n-1)
        if (nCells > 1) {
            var1 = (tempSumsq1[j] - (tempSum1[j] * tempSum1[j]) / nCells) / (nCells - 1);
            var2 = (tempSumsq2[j] - (tempSum2[j] * tempSum2[j]) / nCells) / (nCells - 1);
            var3 = (tempSumsq3[j] - (tempSum3[j] * tempSum3[j]) / nCells) / (nCells - 1);
            // Guard against floating point negatives near zero
            if (var1 < 0) var1 = 0;
            if (var2 < 0) var2 = 0;
            if (var3 < 0) var3 = 0;
            sem1 = sqrt(var1 / nCells);
            sem2 = sqrt(var2 / nCells);
            sem3 = sqrt(var3 / nCells);
        } else {
            sem1 = 0; sem2 = 0; sem3 = 0;
        }

        output += d2s(tempPos[j], 6) + "," +
                  d2s(tempCentered[j], 6) + "," +
                  d2s(mean1, 6) + "," + d2s(mean2, 6) + "," + d2s(mean3, 6) + "," +
                  d2s(sem1, 6) + "," + d2s(sem2, 6) + "," + d2s(sem3, 6) + "," +
                  nCells + "\n";
    }

    // --- Save output ---
    File.saveString(output, outFile);

    showMessage("Done",
        "Grand mean and SEM calculated from normalized per-cell average profiles.\n" +
        "Number of cells used: " + nCells + "\n\n" +
        "Saved to:\n" + outFile);
}


// ============================================================
// Helper function
// ============================================================

// Append a single value to an array (IJ macro arrays are fixed-size, so rebuild)
function appendValue(arr, value) {
    n = lengthOf(arr);
    out = newArray(n + 1);
    for (k = 0; k < n; k++) out[k] = arr[k];
    out[n] = value;
    return out;
}
