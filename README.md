# BACH1landscape

scripts for manuscript on BACH1 invasion landscape

Downoald all files into the same folder.

Unzip the folder data231 so that it is one level below the current folder.

Run the scripts ouBLindLandscUPD.m and ouBHindLandscUPD.m to compare theory, simulations and experimental data for Dox-dependent invasiveness of NF-BL and NF-BH MB231 cells.

Run the scripts gaussSepBL.m and gaussSepBH.m to estimate invasion landscapes NF-BL and NF-BH MB231 cells.

MB231_1.1BLdataGB.xlsx is the BL single-cell fluorescence intensity data, 3 replicates. Columns 1-3: seeded cells. Columns 4-6: invading cells. Each sheet is a different Dox concentration.

MB231_1.1BHdataGB.xlsx is the BH single-cell fluorescence intensity data, 3 replicates. Columns 1-3: seeded cells. Columns 4-6: invading cells. Each sheet is a different Dox concentration.

LandscIndBL.txt is an individual invasion landscape for each lg(BACH1) level.

smLandscBL.txt is a consensus landscape estimated from individual landscapes.

xCentersBL.txt is the set of lg(BACH1) levels for the landscapes.
