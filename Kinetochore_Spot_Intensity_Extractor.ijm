// ============================================================
// Macro: Kinetochore spot intensity extractor
// ============================================================
// Description:
//   Measures a protein-of-interest's intensity at individual
//   kinetochore dots, using an independent poles+kinetochores
//   reference channel to find the dot POSITIONS (so spot detection
//   is not biased toward wherever the protein-of-interest happens
//   to be bright).
//
//   Workflow:
//     1. You assign which channel is the protein-of-interest and
//        which is the poles+kinetochores reference (any other
//        channels present are ignored).
//     2. Find Maxima runs on the REFERENCE channel to detect
//        candidate kinetochore positions.
//     3. Since that reference channel also contains the two
//        spindle poles, the candidate points are added to the ROI
//        Manager so you can quickly delete the pole points by eye
//        before continuing (a few seconds per cell).
//     4. A small circular ROI is placed at each remaining point;
//        mean intensity is measured there in BOTH channels.
//     5. You draw one background ROI for the whole cell (same as
//        the line-profile macros); it's subtracted from every spot.
//     6. A protein/reference ratio is computed per spot - this
//        controls for kinetochore-to-kinetochore and cell-to-cell
//        variability in how well the reference marker itself
//        stained, standard practice for kinetochore intensity work.
//
//   Output is per-SPOT (one row per kinetochore in this cell). Use
//   KT_Spot_PerCell_Summary.ijm afterwards to average all spots in
//   each cell into one per-cell value - kinetochore dots within one
//   cell are not independent replicates, so the cell (not the spot)
//   is the correct unit for statistics, same principle used
//   throughout this pipeline.
//
// Usage:
//   1. Open a multi-channel sum projection image in Fiji/ImageJ.
//   2. Run via Plugins > Macros > Run.
//   3. Assign channels in the dialog.
//   4. Set the Find Maxima prominence and spot radius if the
//      defaults don't suit your staining/magnification.
//   5. Delete any pole points from the ROI Manager when prompted,
//      then click OK to continue.
//   6. Draw a background ROI outside the spindle, click OK.
//   7. Choose an output folder. A CSV is saved automatically.
//
// Output:
//   - <imagename>_KTspots.csv, one row per detected kinetochore:
//     Image, SpotIndex, X_px, Y_px,
//     Protein_raw, Protein_bgsub, Ref_raw, Ref_bgsub, Ratio_bgsub
//
// Notes:
//   - Ratio_bgsub = Protein_bgsub / Ref_bgsub. If Ref_bgsub is at
//     or near zero for a spot (rare, but possible with a very dim
//     kinetochore after background subtraction), the ratio is left
//     blank rather than reporting a distorted or infinite value.
//   - Re-run Find Maxima with a different prominence value (rerun
//     the whole macro) if too many/few candidate points appear -
//     there's no need to fix it by hand beyond removing poles.
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Kinetochore spot intensity extractor" {

    if (nImages == 0) exit("Open an image first.");

    title = getTitle();

    Stack.getDimensions(width, height, channels, slices, frames);
    if (channels < 2) exit("Need at least 2 channels (protein of interest, poles+kinetochores reference).");

    channelList = newArray(channels);
    for (c = 1; c <= channels; c++) channelList[c-1] = "Channel " + c;

    Dialog.create("Assign channels - " + title);
    Dialog.addMessage("Confirm which channel is which for THIS image:");
    Dialog.addChoice("Protein of interest is:", channelList, channelList[0]);
    Dialog.addChoice("Poles + kinetochores reference is:", channelList, channelList[minOf(1, channels-1)]);
    Dialog.addNumber("Find Maxima prominence:", 500);
    Dialog.addNumber("Spot radius (pixels):", 3);
    Dialog.show();

    protChoice = Dialog.getChoice();
    refChoice  = Dialog.getChoice();
    prominence = Dialog.getNumber();
    spotRadius = Dialog.getNumber();

    protCh = parseInt(substring(protChoice, 8));
    refCh  = parseInt(substring(refChoice, 8));

    if (protCh == refCh)
        exit("Protein-of-interest and reference channel can't be the same channel.");

    print("--------------------------------------------------");
    print("Image: " + title);
    print("Assigned: Protein = channel " + protCh + ", Reference = channel " + refCh);
    print("Find Maxima prominence = " + prominence + ", spot radius = " + spotRadius + " px");

    Stack.getPosition(c0, z0, t0);

    // ---- Find Maxima on the reference channel ----
    Stack.setPosition(refCh, z0, 1);
    run("Find Maxima...", "prominence=" + prominence + " output=[Point Selection]");

    if (selectionType() != 10) {
        exit("Find Maxima did not return a point selection - try a different prominence value.");
    }

    getSelectionCoordinates(xpts, ypts);
    nPoints = xpts.length;

    if (nPoints == 0)
        exit("No maxima found - try lowering the prominence value.");

    // ---- add each candidate point to the ROI Manager for manual curation ----
    if (roiManager("count") > 0) roiManager("reset");

    for (i = 0; i < nPoints; i++) {
        makeOval(xpts[i]-2, ypts[i]-2, 4, 4);
        roiManager("Add");
    }

    roiManager("Show All");
    waitForUser("Remove pole points",
        "Reference channel is showing with " + nPoints + " candidate points in the ROI Manager.\n" +
        "Select and DELETE the pole point(s) (usually the 1-2 outermost/brightest) in the\n" +
        "ROI Manager list, then click OK to continue.\n\n" +
        "Leave only genuine kinetochore points selected.");

    nSpots = roiManager("count");
    if (nSpots == 0)
        exit("No kinetochore points left after curation.");

    // ---- background ROI (once per cell) ----
    run("Select None");
    waitForUser("Draw a background ROI outside the spindle, then click OK.");
    if (selectionType() != 0)
        exit("Draw a rectangular background ROI.");

    Stack.setPosition(protCh, z0, 1);
    getStatistics(bgArea1, bgProt);

    Stack.setPosition(refCh, z0, 1);
    getStatistics(bgArea2, bgRef);

    outDir = getDirectory("Choose folder to save CSV files");
    if (outDir == "") exit("No output folder selected.");

    fileBase = replace(title, ".ome", "");
    fileBase = replace(fileBase, ".tif", "");
    fileBase = replace(fileBase, ".tiff", "");

    out = "Image,SpotIndex,X_px,Y_px,Protein_raw,Protein_bgsub,Ref_raw,Ref_bgsub,Ratio_bgsub\n";

    for (s = 0; s < nSpots; s++) {

        roiManager("Select", s);
        Roi.getBounds(bx, by, bw, bh);
        cx = bx + bw/2;
        cy = by + bh/2;

        makeOval(cx - spotRadius, cy - spotRadius, spotRadius*2, spotRadius*2);

        Stack.setPosition(protCh, z0, 1);
        getStatistics(area, protMean);

        Stack.setPosition(refCh, z0, 1);
        getStatistics(area, refMean);

        protBgsub = protMean - bgProt;
        refBgsub  = refMean  - bgRef;

        ratioStr = "";
        if (refBgsub > 0.0001 || refBgsub < -0.0001) {
            ratio = protBgsub / refBgsub;
            ratioStr = "" + ratio;
        }

        out += fileBase + "," + s + "," + cx + "," + cy + "," +
               protMean + "," + protBgsub + "," + refMean + "," + refBgsub + "," + ratioStr + "\n";
    }

    run("Select None");
    roiManager("reset");

    File.saveString(out, outDir + fileBase + "_KTspots.csv");
    print("Saved: " + fileBase + "_KTspots.csv (" + nSpots + " spots)");
    showMessage("Done", nSpots + " kinetochore spots measured.\nSaved to:\n" + outDir + fileBase + "_KTspots.csv");
}
