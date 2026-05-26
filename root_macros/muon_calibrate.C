// muon_calibrate.C
// Compute per-cell MIP (mean energy deposit) from a muon run and produce relative calibration factors.
// Output CSV: /tmp/muon_mip_calibration.csv

#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TMath.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

void muon_calibrate(const char* muonFile = "/home/rudradeb/sim/output/output_thread-1.root",
                    double tmin = 0.0, double tmax = 100.0,
                    bool requireTrackID1 = false,
                    int minHitsPerCell = 10)
{
    const int nSide = 7;
    // allow user-editable file & time window above
    printf("muon_calibrate: file=%s  time window = [%.3f, %.3f]  requireTrackID1=%d\n", muonFile, tmin, tmax, (int)requireTrackID1);

    TFile f(muonFile);
    if (f.IsZombie()) {
        printf("ERROR: cannot open %s\n", muonFile);
        return;
    }
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) {
        printf("ERROR: tree 'Hits' not found in %s\n", muonFile);
        f.Close();
        return;
    }

    // discover time branch name: prefer "fTime" then "time"
    TBranch *btime = (TBranch*)t->GetListOfBranches()->FindObject("fTime");
    const char *timeBranchName = nullptr;
    if (btime) timeBranchName = "fTime";
    else {
        btime = (TBranch*)t->GetListOfBranches()->FindObject("time");
        if (btime) timeBranchName = "time";
    }
    if (!timeBranchName) {
        printf("Warning: no branch named 'fTime' or 'time' found. Setting time to 0 for all hits (no timing cut applied).\n");
    } else {
        printf("Using time branch: '%s'\n", timeBranchName);
    }

    // Branches we expect
    Int_t detID = 0;
    Double_t edep = 0;
    Int_t trackID = 0;
    Double_t timeVal = 0;

    // Set branch addresses if branches exist
    if (t->GetListOfBranches()->FindObject("fDetectorID")) t->SetBranchAddress("fDetectorID", &detID);
    else if (t->GetListOfBranches()->FindObject("detID")) t->SetBranchAddress("detID", &detID);
    else { printf("ERROR: no detector ID branch found (tried 'fDetectorID' and 'detID').\n"); f.Close(); return; }

    if (t->GetListOfBranches()->FindObject("fEdep")) t->SetBranchAddress("fEdep", &edep);
    else if (t->GetListOfBranches()->FindObject("Edep")) t->SetBranchAddress("Edep", &edep);
    else { printf("ERROR: no energy branch found (tried 'fEdep' and 'Edep').\n"); f.Close(); return; }

    if (t->GetListOfBranches()->FindObject("fTrackID")) t->SetBranchAddress("fTrackID", &trackID);
    else if (t->GetListOfBranches()->FindObject("trackID")) t->SetBranchAddress("trackID", &trackID);
    else {
        // not fatal; trackID will remain 0 and user can set requireTrackID1=false
        requireTrackID1 = false;
    }

    if (timeBranchName) t->SetBranchAddress(timeBranchName, &timeVal);
    else timeVal = 0.0;

    // accumulators: sums and counts per (col,row) - combine across layers
    std::vector<std::vector<double>> sumE(nSide, std::vector<double>(nSide, 0.0));
    std::vector<std::vector<long long>> nHits(nSide, std::vector<long long>(nSide, 0));

    Long64_t N = t->GetEntries();
    printf("Total entries in Hits tree: %lld\n", N);

    for (Long64_t i=0; i<N; ++i) {
        t->GetEntry(i);
        if (edep <= 0) continue;
        // timing
        if (timeBranchName) {
            if (timeVal < tmin || timeVal > tmax) continue;
        }
        if (requireTrackID1 && trackID != 1) continue;
        int L, col, row;
        decodeDet(detID, L, col, row);
        if (col < 0 || col >= nSide || row < 0 || row >= nSide) continue;
        sumE[col][row] += edep;
        nHits[col][row] += 1;
    }

    // compute mean per cell
    std::vector<std::vector<double>> meanE(nSide, std::vector<double>(nSide, 0.0));
    double globalSum = 0.0;
    long long globalCount = 0;
    int liveCells = 0;
    for (int x=0;x<nSide;++x) {
        for (int y=0;y<nSide;++y) {
            if (nHits[x][y] >= minHitsPerCell) {
                meanE[x][y] = sumE[x][y] / double(nHits[x][y]);
                globalSum += meanE[x][y];
                globalCount += 1;
                liveCells++;
            } else {
                meanE[x][y] = 0.0; // insufficient stats
            }
        }
    }

    if (globalCount == 0) {
        printf("No cell has >= %d hits. Try reducing minHitsPerCell or check your muon file / timing.\n", minHitsPerCell);
        // nevertheless write CSV with zeros
    }

    double overallMean = (globalCount>0) ? (globalSum / double(globalCount)) : 0.0;
    printf("Cells with >=%d hits: %d  overall mean MIP = %.6e (units of edep branch)\n", minHitsPerCell, liveCells, overallMean);

    // produce CSV: col,row,meanE,calibFactor,countHits
    const char *outfn = "/tmp/muon_mip_calibration.csv";
    FILE *of = fopen(outfn, "w");
    if (!of) { printf("ERROR: cannot open %s for write\n", outfn); f.Close(); return; }
    fprintf(of, "col,row,meanE,calibFactor,countHits\n");
    for (int y=0;y<nSide;++y) {
        for (int x=0;x<nSide;++x) {
            double m = meanE[x][y];
            double factor = (m > 0 && overallMean > 0) ? (overallMean / m) : 0.0;
            long long hits = nHits[x][y];
            fprintf(of, "%d,%d,%.12e,%.8f,%lld\n", x, y, m, factor, hits);
        }
    }
    fclose(of);
    printf("Wrote per-cell calibration file: %s\n", outfn);

    // print a nice table to screen
    printf("\nPer-cell MIP means (rows y=0..6, cols x=0..6):\n");
    for (int y=0;y<nSide;++y) {
        for (int x=0;x<nSide;++x) {
            if (nHits[x][y] >= minHitsPerCell) printf("%12.6e ", meanE[x][y]);
            else printf("   (nan)      ");
        }
        printf("\n");
    }
    printf("\nCalibration factor = overallMean / meanCell (1.0 = no correction). Overall mean = %.6e\n", overallMean);

    f.Close();
    printf("muon_calibrate finished.\n");
}

