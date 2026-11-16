// ============================================================
// Macro: Kinetochore spot per-cell summary
// ============================================================
// Description:
//   Reads a folder of *_KTspots.csv files (one per cell, produced
//   by Kinetochore_Spot_Intensity_Extractor.ijm) and averages all
//   kinetochore spots WITHIN each cell into a single row per cell.
//   This is the correct unit for statistics - individual
//   kinetochore dots within one cell are not independent
//   replicates, same principle used throughout this pipeline.
//
// Usage:
//   1. Place all *_KTspots.csv files (both conditions, all
//      replicates) into a single folder.
//   2. Run via Plugins > Macros > Run.
//   3. Select that folder.
//   4. Enter the two condition keywords when prompted (whatever
//      substrings identify each condition in your filenames -
//      e.g. "dmso" and "importazole", or anything else).
//
// Input:
//   - Folder of *_KTspots.csv files with columns:
//     Image, SpotIndex, X_px, Y_px,
//     Protein_raw, Protein_bgsub, Ref_raw, Ref_bgsub, Ratio_bgsub
//
// Output (saved in same folder):
//   - percell_KTspot_summary.csv
//   - Same table also pushed to the Results window / .xls
//
// Output CSV columns (one row per cell, sorted by condition):
//   Filename, Condition, N_kinetochores,
//   MeanProtein_bgsub, MeanRef_bgsub, MeanRatio_bgsub
//
// Notes:
//   - MeanRatio_bgsub averages only the spots that had a valid
//     ratio (i.e. excludes any spot where Ref_bgsub was too close
//     to zero to divide by, left blank by the extraction macro).
//   - Files whose name doesn't contain either condition keyword
//     are skipped, with a note in the Log.
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Kinetochore spot per-cell summary" {

    dir = getDirectory("Choose the folder containing *_KTspots.csv files");
    if (dir=="") exit("No folder selected.");

    Dialog.create("Condition keywords");
    Dialog.addMessage("For each condition, you can enter more than one filename\nkeyword separated by commas (e.g. \"mis12, sim12\") if your\nfiles use more than one spelling.");
    Dialog.addString("Condition A name:", "ConditionA");
    Dialog.addString("Condition A filename keyword(s):", "conditionA");
    Dialog.addString("Condition B name:", "ConditionB");
    Dialog.addString("Condition B filename keyword(s):", "conditionB");
    Dialog.show();

    condAName = Dialog.getString();
    condAKeysRaw = toLowerCase(Dialog.getString());
    condBName = Dialog.getString();
    condBKeysRaw = toLowerCase(Dialog.getString());

    condAKeys = split(condAKeysRaw, ",");
    condBKeys = split(condBKeysRaw, ",");
    for (i = 0; i < condAKeys.length; i++) condAKeys[i] = String.trim(condAKeys[i]);
    for (i = 0; i < condBKeys.length; i++) condBKeys[i] = String.trim(condBKeys[i]);

    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }

    print("\\Clear");
    print("Kinetochore spot per-cell summary for: " + dir);
    print("Condition A: " + condAName + " <- " + condAKeysRaw);
    print("Condition B: " + condBName + " <- " + condBKeysRaw);
    print("--------------------------------------------------");

    list = getFileList(dir);
    maxRows = list.length;

    cFilename = newArray(maxRows);
    cCond     = newArray(maxRows);
    cN        = newArray(maxRows);
    cMeanProt = newArray(maxRows);
    cMeanRef  = newArray(maxRows);
    cMeanRatio = newArray(maxRows);
    nCollected = 0;

    for (f = 0; f < list.length; f++) {

        name = list[f];
        if (endsWith(name, "/") || endsWith(name, "\\")) continue;

        lower = toLowerCase(name);

        if (!endsWith(lower, "_ktspots.csv")) continue;
        if (indexOf(lower, "percell") != -1) continue;

        cond = "";
        for (ki = 0; ki < condAKeys.length; ki++) {
            if (lengthOf(condAKeys[ki]) > 0 && indexOf(lower, condAKeys[ki]) != -1) { cond = condAName; ki = condAKeys.length; }
        }
        if (cond == "") {
            for (ki = 0; ki < condBKeys.length; ki++) {
                if (lengthOf(condBKeys[ki]) > 0 && indexOf(lower, condBKeys[ki]) != -1) { cond = condBName; ki = condBKeys.length; }
            }
        }
        if (cond == "") {
            print("Skipped (no condition keyword matched): " + name);
            continue;
        }

        path = dir + name;
        text = File.openAsString(path);
        if (lengthOf(text) == 0) {
            print("Skipped empty/unreadable file: " + path);
            continue;
        }

        lines = split(text, "\n");

        sumProt = 0; nProt = 0;
        sumRef  = 0; nRef  = 0;
        sumRatio = 0; nRatio = 0;

        for (i = 1; i < lines.length; i++) {

            line = replace(lines[i], "\r", "");
            if (lengthOf(line) == 0) continue;

            cols = split(line, ",");
            if (cols.length < 9) continue;

            protBgsub = parseFloat(cols[5]);
            refBgsub  = parseFloat(cols[7]);

            if (!isNaN(protBgsub)) { sumProt += protBgsub; nProt++; }
            if (!isNaN(refBgsub))  { sumRef  += refBgsub;  nRef++;  }

            if (cols.length >= 9 && lengthOf(cols[8]) > 0) {
                ratioVal = parseFloat(cols[8]);
                if (!isNaN(ratioVal)) { sumRatio += ratioVal; nRatio++; }
            }
        }

        if (nProt == 0) {
            print("Skipped (no valid spot rows): " + name);
            continue;
        }

        cFilename[nCollected] = name;
        cCond[nCollected]     = cond;
        cN[nCollected]        = nProt;
        cMeanProt[nCollected] = sumProt / nProt;
        cMeanRef[nCollected]  = sumRef / nRef;
        if (nRatio > 0) {
            cMeanRatio[nCollected] = sumRatio / nRatio;
        } else {
            cMeanRatio[nCollected] = NaN;
        }
        nCollected++;

        print("Used: " + name + "   [" + cond + "]   n_kinetochores=" + nProt);
    }

    if (nCollected == 0)
        exit("No matching files found.");

    // ---- sort: Condition A first, then Condition B, then filename ----
    condPriority = newArray(condAName, condBName);

    sortKey = newArray(nCollected);
    for (r = 0; r < nCollected; r++) {
        p = -1;
        for (pi = 0; pi < condPriority.length; pi++) {
            if (condPriority[pi] == cCond[r]) { p = pi; pi = condPriority.length; }
        }
        if (p == -1) p = 99;
        padded = "" + p;
        if (p < 10) padded = "0" + padded;
        sortKey[r] = padded + "_" + cFilename[r];
    }

    ranks = Array.rankPositions(sortKey);

    out = "Filename,Condition,N_kinetochores,MeanProtein_bgsub,MeanRef_bgsub,MeanRatio_bgsub\n";

    for (k = 0; k < nCollected; k++) {
        r = ranks[k];

        ratioOut = "" + cMeanRatio[r];
        if (isNaN(cMeanRatio[r])) ratioOut = "";

        out += cFilename[r] + "," + cCond[r] + "," + cN[r] + "," +
               cMeanProt[r] + "," + cMeanRef[r] + "," + ratioOut + "\n";

        setResult("Filename", k, cFilename[r]);
        setResult("Condition", k, cCond[r]);
        setResult("N_kinetochores", k, cN[r]);
        setResult("MeanProtein_bgsub", k, cMeanProt[r]);
        setResult("MeanRef_bgsub", k, cMeanRef[r]);
        if (!isNaN(cMeanRatio[r])) setResult("MeanRatio_bgsub", k, cMeanRatio[r]);
    }

    File.saveString(out, dir + "percell_KTspot_summary.csv");
    updateResults();
    saveAs("Results", dir + "percell_KTspot_summary.xls");

    print("--------------------------------------------------");
    print("DONE");
    print("Cells summarized = " + nCollected);
    print("Saved: " + dir + "percell_KTspot_summary.csv");
    print("Saved: " + dir + "percell_KTspot_summary.xls");
}
