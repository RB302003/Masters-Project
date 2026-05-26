#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

// Simple diagnostic macro for your Hits tree.
// Prints per-event peak-layer centroid and energy info and summary stats.
//
// Usage: compile & run in ROOT:
//   root -l diagnoseEvents.C
//   diagnoseEvents();   // defaults: inspect all four energy files listed
//
// If your branch names differ, change the SetBranchAddress names accordingly.

void diagnoseEvents()
{
    // EDIT ONLY if file paths are different
    std::map<int,std::string> fileOf = {
        {1,  "/home/rudradeb/sim/finaloutput/output_thread-1.root"},
        {2,  "/home/rudradeb/sim/finaloutput/output_thread-2.root"},
        {5,  "/home/rudradeb/sim/finaloutput/output_thread-3.root"},
        {10, "/home/rudradeb/sim/finaloutput/output_thread-4.root"}
    };
    std::vector<int> energies = {1,2,5,10};

    const int nLayers = 24;
    const int nSide = 7;
    const int outCenter = 3; // detector center index that's expected

    auto decode = [](int detID, int &layer, int &col, int &row){
        layer = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
        int tmp = detID % 10000;
        col = tmp / 100;
        row = tmp % 100;
    };

    const int Nprint = 20; // how many events to print in detail per file
    for (int E : energies) {
        std::string fname = fileOf[E];
        std::cout << "\n=== DIAGNOSTICS: " << E << " GeV  file: " << fname << " ===\n";

        TFile f(fname.c_str());
        if (f.IsZombie()) { std::cout<<"  cannot open file\n"; continue; }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) { std::cout<<"  no Hits tree\n"; f.Close(); continue; }

        // Branches (match your tree)
        Int_t detID, trackID, evID_branch;
        Double_t edep;
        t->SetBranchAddress("fDetectorID",&detID);
        t->SetBranchAddress("fEdep",&edep);
        t->SetBranchAddress("fTrackID",&trackID);
        t->SetBranchAddress("fEvent",&evID_branch);

        // first pass: accumulate per-event, per-layer 7x7
        std::map<Int_t, std::vector<std::vector<std::vector<double>>>> perEvent; // ev -> [layer][x][y]
        Long64_t nEntries = t->GetEntries();
        for (Long64_t i=0;i<nEntries;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            int L,cx,ry; decode(detID,L,cx,ry);
            if (L<0 || L>=nLayers) continue;
            if (cx<0 || cx>=nSide || ry<0 || ry>=nSide) continue;
            if (perEvent.find(evID_branch) == perEvent.end())
                perEvent[evID_branch] = std::vector<std::vector<std::vector<double>>>(nLayers,
                    std::vector<std::vector<double>>(nSide, std::vector<double>(nSide,0.0)));
            perEvent[evID_branch][L][cx][ry] += edep;
        }
        f.Close();

        // Now inspect per-event summary
        int ie=0;
        double sumCx=0, sumCy=0, sumCx2=0, sumCy2=0;
        Long64_t totalEvents = 0;
        double totalE_allEvents = 0.0;
        double totalE_peakLayerSum = 0.0;
        double totalE_accum_if_centered = 0.0; // estimate of energy kept after centering into 7x7
        std::vector<int> centroidBinCounts(nSide*nSide,0);

        // open file again for track-split (we will need trackID to compute pri/sec if desired)
        TFile f2(fname.c_str());
        if (f2.IsZombie()) { std::cout<<"  cannot re-open file for track splits\n"; continue; }
        TTree *t2 = (TTree*)f2.Get("Hits");
        if (!t2) { std::cout<<"  no Hits tree\n"; f2.Close(); continue; }

        // Map event->vector of hits indices could be heavy; instead for quick diagnostics we'll compute per-event maps
        // of total energy and centroid using perEvent stored above (we only used perEvent).
        for (auto &p : perEvent) {
            Int_t eid = p.first;
            auto &layerCell = p.second;
            // compute per-layer totals
            std::vector<double> layerSum(nLayers,0.0);
            for (int L=0; L<nLayers; ++L)
                for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) layerSum[L] += layerCell[L][x][y];
            int peakL = std::max_element(layerSum.begin(), layerSum.end()) - layerSum.begin();
            double peakTotal = layerSum[peakL];
            if (peakTotal <= 0) continue; // skip empty
            // per-event 7x7 map from peak layer
            std::vector<std::vector<double>> map7(nSide, std::vector<double>(nSide,0.0));
            for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) map7[x][y] = layerCell[peakL][x][y];
            // centroid & max cell
            double sE=0, sX=0, sY=0; int maxX=0, maxY=0; double maxV=-1;
            for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) {
                double v = map7[x][y];
                sE += v; sX += v*x; sY += v*y;
                if (v>maxV) { maxV=v; maxX=x; maxY=y; }
            }
            if (sE<=0) continue;
            double cx = sX/sE;
            double cy = sY/sE;
            // if we center to (3,3) by integer shift, fraction kept = sum of cells after shift inside 7x7 / sE
            int shiftX = (int)TMath::Nint(outCenter - cx);
            int shiftY = (int)TMath::Nint(outCenter - cy);
            double kept=0;
            for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) {
                int tx = x + shiftX;
                int ty = y + shiftY;
                if (tx<0||tx>=nSide||ty<0||ty>=nSide) continue;
                kept += map7[x][y];
            }
            double fracKept = (sE>0) ? (kept / sE) : 0.0;

            // collect stats
            totalEvents++;
            sumCx += cx; sumCy += cy; sumCx2 += cx*cx; sumCy2 += cy*cy;
            totalE_allEvents += sE;
            totalE_peakLayerSum += peakTotal;
            totalE_accum_if_centered += kept;
            // centroid bin (coarse)
            int bx = (int)TMath::FloorNint(cx); if (bx<0) bx=0; if (bx>=nSide) bx=nSide-1;
            int by = (int)TMath::FloorNint(cy); if (by<0) by=0; if (by>=nSide) by=nSide-1;
            centroidBinCounts[by*nSide + bx]++;

            // print first Nprint events details
            if (ie < Nprint) {
                printf("event %6d  peakL=%2d  totE=%.6f  cx=%.3f cy=%.3f  max=(%d,%d)=%.6f  keptFrac=%.3f\n",
                       eid, peakL, sE, cx, cy, maxX, maxY, maxV, fracKept);
            }
            ie++;
        }

        // summary
        double meanCx = (totalEvents>0) ? sumCx/totalEvents : 0.0;
        double meanCy = (totalEvents>0) ? sumCy/totalEvents : 0.0;
        double rmsCx = (totalEvents>0) ? sqrt( sumCx2/totalEvents - meanCx*meanCx ) : 0.0;
        double rmsCy = (totalEvents>0) ? sqrt( sumCy2/totalEvents - meanCy*meanCy ) : 0.0;
        double fracEnergyKept = (totalE_peakLayerSum>0) ? (totalE_accum_if_centered / totalE_peakLayerSum) : 0.0;

        std::cout << "\nSUMMARY " << E << " GeV:\n";
        std::cout << "  events considered = " << totalEvents << "\n";
        std::cout << "  mean centroid (cx,cy) = (" << meanCx << ", " << meanCy << ")\n";
        std::cout << "  rms centroid (cx,cy)  = (" << rmsCx << ", " << rmsCy << ")\n";
        std::cout << "  total energy in peak-layers (sum over events) = " << totalE_peakLayerSum << "\n";
        std::cout << "  total energy that would be kept after integer-centering = " << totalE_accum_if_centered << "\n";
        std::cout << "  fraction kept (kept / peakSum) = " << fracEnergyKept << "\n";

        // centroid occupancy extremes
        int minBin=INT_MAX, maxBin=0, minIdx=-1, maxIdx=-1;
        for (int i=0;i<nSide*nSide;++i) {
            if (centroidBinCounts[i] < minBin) { minBin = centroidBinCounts[i]; minIdx = i; }
            if (centroidBinCounts[i] > maxBin) { maxBin = centroidBinCounts[i]; maxIdx = i; }
        }
        int minbx = minIdx % nSide, minby = minIdx / nSide;
        int maxbx = maxIdx % nSide, maxby = maxIdx / nSide;
        std::cout << "  centroid occupancy: min bin ("<<minbx<<","<<minby<<")="<<minBin
                  << "  max bin ("<<maxbx<<","<<maxby<<")="<<maxBin << "\n";

        f2.Close();

        // write a tiny csv in /tmp for further inspection if you want (eventID,peakLayer,cx,cy,totalE,fracKept)
        std::string outcsv = Form("/tmp/diagnose_events_%dGeV.csv", E);
        FILE *of = fopen(outcsv.c_str(),"w");
        if (of) {
            fprintf(of,"event,peakLayer,cx,cy,totalE,fracKept\n");
            int printed=0;
            for (auto &p : perEvent) {
                Int_t eid = p.first;
                auto &layerCell = p.second;
                std::vector<double> layerSum(nLayers,0.0);
                for (int L=0;L<nLayers;++L) for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) layerSum[L] += layerCell[L][x][y];
                int peakL = std::max_element(layerSum.begin(), layerSum.end()) - layerSum.begin();
                double peakTotal = layerSum[peakL];
                if (peakTotal <= 0) continue;
                std::vector<std::vector<double>> map7(nSide,std::vector<double>(nSide,0.0));
                for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) map7[x][y] = layerCell[peakL][x][y];
                double sE=0,sX=0,sY=0; for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) { double v=map7[x][y]; sE+=v; sX+=v*x; sY+=v*y;}
                if (sE<=0) continue;
                double cx = sX/sE, cy = sY/sE;
                int shiftX = (int)TMath::Nint(outCenter - cx), shiftY = (int)TMath::Nint(outCenter - cy);
                double kept=0; for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y){ int tx = x+shiftX, ty=y+shiftY; if (tx<0||tx>=nSide||ty<0||ty>=nSide) continue; kept += map7[x][y]; }
                double fracKept = (sE>0) ? kept/sE : 0.0;
                fprintf(of,"%d,%d,%.6f,%.6f,%.6f,%.6f\n", eid, peakL, cx, cy, sE, fracKept);
                if (++printed > 10000) break;
            }
            fclose(of);
            std::cout << "  wrote /tmp/diagnose_events_" << E << "GeV.csv (first 10000 events)\n";
        } else {
            std::cout << "  couldn't write csv\n";
        }

        // free perEvent map to release memory
        perEvent.clear();
    } // energies
    std::cout << "\nDIAGNOSTICS DONE\n";
}

