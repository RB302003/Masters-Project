// electron_to_mips_fixed.C
// Corrected: compute global scale using total MIPs across run AND number of events.
// Outputs saved in /home/rudradeb/sim/output/

#include <map>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

// electronFile: electron ROOT file (10 GeV run by default)
// muCSV: muon calibration CSV produced by muon_calibrate
// targetEnergyGeV: beam energy per event (e.g. 10.0)
void electron_to_mips(const char* electronFile = "/home/rudradeb/sim/finaloutput/output_thread-4.root",
                            const char* muCSV = "/tmp/muon_mip_calibration.csv",
                            double targetEnergyGeV = 10.0)
{
    const int nSide = 7;

    // 1) read muon CSV into map (col,row)->meanMIP
    std::map<std::pair<int,int>, double> muMean;
    {
        std::ifstream ifs(muCSV);
        if (!ifs.is_open()) {
            printf("ERROR: cannot open muon CSV: %s\n", muCSV);
            return;
        }
        std::string line;
        if (!std::getline(ifs, line)) { printf("ERROR: empty muon CSV\n"); ifs.close(); return; }
        while (std::getline(ifs, line)) {
            if (line.size() == 0) continue;
            std::stringstream ss(line);
            int col=0, row=0;
            double meanE=0, calib=0; long long cnt=0;
            char comma;
            if (!(ss >> col)) continue;
            ss >> comma;
            if (!(ss >> row)) continue;
            ss >> comma;
            if (!(ss >> meanE)) meanE = 0.0;
            muMean[{col,row}] = meanE;
        }
        ifs.close();
    }
    int found=0;
    for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) if (muMean.find({x,y}) != muMean.end()) ++found;
    printf("Read muon calibration: found %d / %d cells.\n", found, nSide*nSide);

    // 2) open electron file and sum total energy and total MIPs across the run
    TFile f(electronFile);
    if (f.IsZombie()) { printf("ERROR: cannot open electron file: %s\n", electronFile); return; }
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) { printf("ERROR: no Hits tree in %s\n", electronFile); f.Close(); return; }

    // Branch discovery
    Int_t detID = 0;
    Double_t edep = 0;
    Int_t evID = 0;
    bool hasEventBranch = false;
    if (t->GetListOfBranches()->FindObject("fDetectorID")) t->SetBranchAddress("fDetectorID",&detID);
    else if (t->GetListOfBranches()->FindObject("detID")) t->SetBranchAddress("detID",&detID);
    else { printf("ERROR: no detector id branch\n"); f.Close(); return; }
    if (t->GetListOfBranches()->FindObject("fEdep")) t->SetBranchAddress("fEdep",&edep);
    else if (t->GetListOfBranches()->FindObject("Edep")) t->SetBranchAddress("Edep",&edep);
    else { printf("ERROR: no energy branch\n"); f.Close(); return; }

    // event id branch optional - try multiple names
    if (t->GetListOfBranches()->FindObject("fEvent")) { t->SetBranchAddress("fEvent",&evID); hasEventBranch=true; }
    else if (t->GetListOfBranches()->FindObject("event")) { t->SetBranchAddress("event",&evID); hasEventBranch=true; }
    else if (t->GetListOfBranches()->FindObject("evt")) { t->SetBranchAddress("evt",&evID); hasEventBranch=true; }
    else hasEventBranch = false;

    Long64_t N = t->GetEntries();
    printf("Electron file has %lld hit entries. hasEventBranch=%d\n", N, (int)hasEventBranch);

    // accumulators
    std::vector<std::vector<double>> totalE(nSide, std::vector<double>(nSide,0.0));
    double total_MIPs_sum = 0.0;
    double totalE_allcells = 0.0;
    std::set<int> eventIDs;

    for (Long64_t i=0;i<N;++i) {
        t->GetEntry(i);
        if (edep <= 0) continue;
        int L, cx, ry; decodeDet(detID, L, cx, ry);
        if (cx<0 || cx>=nSide || ry<0 || ry>=nSide) continue;
        totalE[cx][ry] += edep;
        totalE_allcells += edep;
        double meanMIP = 0.0;
        auto it = muMean.find({cx,ry});
        if (it != muMean.end()) meanMIP = it->second;
        // If the muon mean for this cell is zero (missing), use the average of muMean map:
        if (meanMIP <= 0.0) {
            // compute fallback average on demand (cheap since map small)
            double s=0; int c=0;
            for (auto &kv : muMean) if (kv.second>0) { s += kv.second; ++c; }
            meanMIP = (c>0) ? (s/double(c)) : 1.0;
        }
        double nMips_hit = (meanMIP>0) ? (edep / meanMIP) : 0.0;
        total_MIPs_sum += nMips_hit;

        if (hasEventBranch) eventIDs.insert((int)evID);
    }

    // if no event branch, try to estimate N_events from other place:
    size_t N_events = hasEventBranch ? eventIDs.size() : 0;
    if (!hasEventBranch) {
        // Attempt a heuristic: if there is a branch 'fEvent' missing, but there is a top-level tree with event count,
        // user can pass number of events manually in the function arguments; for now warn and assume 1 event (not ideal).
        printf("WARNING: no event id branch found. Cannot determine N_events automatically; assuming N_events=1 (you should rerun with event-branch available).\n");
        N_events = 1;
    }

    printf("Run totals: totalE_allcells=%.6e (units of fEdep), total_MIPs_sum=%.12e, N_events=%zu\n",
           totalE_allcells, total_MIPs_sum, N_events);

    if (N_events == 0) { printf("ERROR: N_events==0 — aborting\n"); f.Close(); return; }

    // compute correct scale: targetEnergyGeV * N_events / total_MIPs_sum
    double scale_GeV_per_MIP = (total_MIPs_sum > 0.0) ? ( (targetEnergyGeV * double(N_events)) / total_MIPs_sum ) : 0.0;

    // write outputs into /home/rudradeb/sim/output/
    const char *outcsv = "/home/rudradeb/sim/output/electron_mips_and_scale.csv";
    const char *outtxt = "/home/rudradeb/sim/output/electron_scale.txt";

    FILE *of = fopen(outcsv,"w");
    if (!of) { printf("ERROR: cannot open %s for write\n", outcsv); f.Close(); return; }
    fprintf(of,"col,row,totalE,N_MIPs\n");
    // Recompute per-cell N_MIPs using totals computed above
    for (int y=0;y<nSide;++y) {
        for (int x=0;x<nSide;++x) {
            double meanMIP = 0.0;
            auto it = muMean.find({x,y});
            if (it != muMean.end()) meanMIP = it->second;
            if (meanMIP <= 0.0) {
                double s=0; int c=0;
                for (auto &kv : muMean) if (kv.second>0) { s+=kv.second; ++c; }
                meanMIP = (c>0)?(s/double(c)):1.0;
            }
            double nm = (meanMIP>0) ? (totalE[x][y] / meanMIP) : 0.0;
            fprintf(of, "%d,%d,%.12e,%.12e\n", x, y, totalE[x][y], nm);
        }
    }
    fclose(of);
    printf("Wrote per-cell CSV: %s\n", outcsv);

    FILE *ot = fopen(outtxt, "w");
    if (ot) {
        fprintf(ot, "targetEnergyGeV, N_events, total_MIPs_sum, scale_GeV_per_MIP\n");
        fprintf(ot, "%.6f, %zu, %.12e, %.12e\n", targetEnergyGeV, N_events, total_MIPs_sum, scale_GeV_per_MIP);
        fclose(ot);
        printf("Wrote scale file: %s\n", outtxt);
    } else {
        printf("Warning: could not write %s\n", outtxt);
    }

    printf("Scale computed: 1 MIP = %.12e GeV  (using targetEnergy=%.3f GeV and N_events=%zu)\n",
           scale_GeV_per_MIP, targetEnergyGeV, N_events);

    f.Close();
    printf("electron_to_mips_fixed finished.\n");
}

