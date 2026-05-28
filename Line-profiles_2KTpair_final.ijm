// ============================================================
// Macro: Kinetochore-to-Kinetochore Line Profile Extractor
// ============================================================
// Description:
//   Extracts fluorescence intensity line profiles along the
//   kinetochore-to-kinetochore (KT-KT) axis for two kinetochore
//   pairs per cell, across three channels simultaneously.
//
//   For each cell, the user draws a straight line between two
//   sister kinetochores for two representative KT pairs. Profiles
//   are resampled to a common number of bins, min-max normalized,
//   and the per-cell average of the two pairs is also computed.
//
//   Output is a CSV file containing:
//     - Raw and normalized intensity profiles for each pair
//     - Per-cell averaged profile (pair_type = "average")
//     - Position along the KT-KT axis (0 to 1, and centered -0.5 to 0.5)
//
// Usage:
//   1. Open a multi-channel z-stack in Fiji/ImageJ.
//   2. Run this macro via Plugins > Macros > Run.
//   3. Enter channel numbers for kinetochore marker and two proteins.
//   4. For each cell and each KT pair:
//      - Navigate to the correct z-slice
//      - Draw a straight line between the two kinetochores
//      - Click OK in the dialog
//   5. Choose an output folder. A CSV file will be saved automatically.
//
// Input:
//   - A multi-channel z-stack (e.g., .tif, .tiff, or .dv format)
//   - Channel 1: kinetochore marker (e.g., CENP-A)
//   - Channel 2: protein of interest 1 (e.g., HURP)
//   - Channel 3: protein of interest 2 (e.g., PRC1/tubulin)
//
// Output:
//   - <imagename>_profiles.csv with columns:
//     image, cell, pair_type, position_0to1, centered_pos,
//     ch1_raw, ch2_raw, ch3_raw, ch1_norm, ch2_norm, ch3_norm
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// Publication: Skendo et al., 2025 (Cell Reports, title to be confirmed)
// Repository: https://github.com/kskendo/KT-LineProfile-Macros
// ============================================================

macro "Kinetochore normalized profiles per cell FAST" {

    if (nImages == 0) exit("Open an image first.");

    title = getTitle();
    Stack.getDimensions(w, h, channels, slices, frames);

    // --- User inputs ---
    ch1 = getNumber("Channel for kinetochore signal (e.g., CENP-A):", 1);
    ch2 = getNumber("Channel for protein 1:", 2);
    ch3 = getNumber("Channel for protein 2:", 3);

    nBins = getNumber("Number of normalized distance bins:", 101);
    nCells = getNumber("How many cells do you want to analyze in this image?", 1);
    lineWidth = getNumber("Line width (pixels):", 1);

    outDir = getDirectory("Choose output folder");
    if (outDir == "") exit("No output folder selected.");

    run("Line Width...", "line=" + lineWidth);

    // --- Prepare base filename ---
    baseName = title;
    baseName = replace(baseName, ".tif", "");
    baseName = replace(baseName, ".tiff", "");
    baseName = replace(baseName, ".dv", "");
    baseName = replace(baseName, " ", "_");

    outFile = outDir + baseName + "_profiles.csv";

    // --- Build normalized position arrays ---
    pos = newArray(nBins);
    centered = newArray(nBins);
    for (i = 0; i < nBins; i++) {
        pos[i] = i / (nBins - 1);        // 0 to 1
        centered[i] = pos[i] - 0.5;      // -0.5 to 0.5 (centered at KT midpoint)
    }

    // --- Build output CSV in memory ---
    csv = "image,cell,pair_type,position_0to1,centered_pos,"
        + "ch1_raw,ch2_raw,ch3_raw,"
        + "ch1_norm,ch2_norm,ch3_norm\n";

    // --- Main loop: iterate over cells ---
    for (cell = 1; cell <= nCells; cell++) {

        // Arrays to store profiles for both KT pairs
        p1c1 = newArray(nBins); p1c2 = newArray(nBins); p1c3 = newArray(nBins);
        p2c1 = newArray(nBins); p2c2 = newArray(nBins); p2c3 = newArray(nBins);

        for (pair = 1; pair <= 2; pair++) {

            // Prompt user to draw line between kinetochores
            waitForUser("Cell " + cell + " — Pair " + pair,
                "Navigate to the correct Z slice.\n" +
                "Draw a straight line between the 2 kinetochores.\n" +
                "Then click OK.");

            if (selectionType() != 5)
                exit("Please draw a straight line selection.");

            // Save current stack position
            Stack.getPosition(c0, z0, t0);

            // Extract profiles from each channel (batch mode for speed)
            setBatchMode(true);

            Stack.setPosition(ch1, z0, t0);
            prof1 = getProfile();
            prof1R = resampleArray(prof1, nBins);

            Stack.setPosition(ch2, z0, t0);
            prof2 = getProfile();
            prof2R = resampleArray(prof2, nBins);

            Stack.setPosition(ch3, z0, t0);
            prof3 = getProfile();
            prof3R = resampleArray(prof3, nBins);

            // Restore original channel position
            Stack.setPosition(c0, z0, t0);
            setBatchMode(false);

            // Store raw profiles per pair
            if (pair == 1) {
                copyArray(prof1R, p1c1); copyArray(prof2R, p1c2); copyArray(prof3R, p1c3);
            } else {
                copyArray(prof1R, p2c1); copyArray(prof2R, p2c2); copyArray(prof3R, p2c3);
            }

            // Min-max normalize each profile
            norm1 = minMaxNormalize(prof1R);
            norm2 = minMaxNormalize(prof2R);
            norm3 = minMaxNormalize(prof3R);

            // Append pair rows to CSV
            csv = csv + rowsToString(baseName, cell, "pair" + pair, pos, centered,
                prof1R, prof2R, prof3R, norm1, norm2, norm3);
        }

        // Compute per-cell average of the two KT pairs
        avg1 = averageTwoArrays(p1c1, p2c1);
        avg2 = averageTwoArrays(p1c2, p2c2);
        avg3 = averageTwoArrays(p1c3, p2c3);

        avgNorm1 = minMaxNormalize(avg1);
        avgNorm2 = minMaxNormalize(avg2);
        avgNorm3 = minMaxNormalize(avg3);

        // Append per-cell average row
        csv = csv + rowsToString(baseName, cell, "average", pos, centered,
            avg1, avg2, avg3, avgNorm1, avgNorm2, avgNorm3);
    }

    // --- Save output ---
    File.saveString(csv, outFile);
    showMessage("Done", "Profiles exported to:\n" + outFile);
}


// ============================================================
// Helper functions
// ============================================================

// Resample array 'a' to 'nOut' evenly spaced points using linear interpolation
function resampleArray(a, nOut) {
    nIn = lengthOf(a);
    out = newArray(nOut);
    if (nIn == 1) {
        for (i = 0; i < nOut; i++) out[i] = a[0];
        return out;
    }
    for (i = 0; i < nOut; i++) {
        x = i * (nIn - 1) / (nOut - 1);
        i0 = floor(x);
        i1 = i0 + 1;
        if (i1 >= nIn) i1 = nIn - 1;
        frac = x - i0;
        out[i] = a[i0] * (1 - frac) + a[i1] * frac;
    }
    return out;
}

// Average two arrays element-wise
function averageTwoArrays(a, b) {
    n = lengthOf(a);
    out = newArray(n);
    for (i = 0; i < n; i++) out[i] = (a[i] + b[i]) / 2;
    return out;
}

// Min-max normalize array to [0, 1]
function minMaxNormalize(a) {
    n = lengthOf(a);
    out = newArray(n);
    minV = a[0]; maxV = a[0];
    for (i = 1; i < n; i++) {
        if (a[i] < minV) minV = a[i];
        if (a[i] > maxV) maxV = a[i];
    }
    range = maxV - minV;
    if (range == 0) {
        for (i = 0; i < n; i++) out[i] = 0;
        return out;
    }
    for (i = 0; i < n; i++) out[i] = (a[i] - minV) / range;
    return out;
}

// Copy source array into destination array
function copyArray(src, dst) {
    n = lengthOf(src);
    for (i = 0; i < n; i++) dst[i] = src[i];
}

// Build CSV rows for one profile entry
function rowsToString(imageName, cellID, pairType, pos, centered,
    ch1raw, ch2raw, ch3raw, ch1norm, ch2norm, ch3norm) {
    n = lengthOf(pos);
    s = "";
    for (i = 0; i < n; i++) {
        s = s + imageName + "," + cellID + "," + pairType + "," +
            d2s(pos[i], 6) + "," + d2s(centered[i], 6) + "," +
            d2s(ch1raw[i], 6) + "," + d2s(ch2raw[i], 6) + "," + d2s(ch3raw[i], 6) + "," +
            d2s(ch1norm[i], 6) + "," + d2s(ch2norm[i], 6) + "," + d2s(ch3norm[i], 6) + "\n";
    }
    return s;
}
