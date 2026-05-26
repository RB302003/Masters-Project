// exportFinalMaps.C
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH2D.h"

// Small utility to decode your packed detID: layer*10000 + col*100 + row
static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

// exportFinalMaps(centerEvents=true, makeFlat=true, perEventNormalize=true)
// - computes final 7x7 averaged arrays and prints & saves CSV for Tot/Pri/Sec
void exportFinalMaps(Bool_t centerEvents = kTRUE, Bool_t makeFlat = kTRUE, Bool_t perEventNormalize = kTRUE)
{
    // --- file list: edit if needed ---
    std::map<int,std::string> fileOf = {
        {1,  "/home/rudradeb/sim/finaloutput/output_thread-1.root"},
        {2,  "/home/rudradeb/sim/finaloutput/output_thread-2.root"},
        {5,  "/home/rudradeb/sim/finaloutput/output_thread-3.root"},
        {10, "/home/rudradeb/sim/finaloutput/output_thread-4.root"}
    };
    std::vector<int> energies = {1,2,5,10};

    const int nLayers = 24;
    const int outSide = 7;
    const int eventSide = 11;
    const int eventCenter = (eventSide-1)/2;
    const int halfOut = (outSide-1)/2;

    // flattening params (safe defaults)
    const int centroidBinsPerSide = 12;
    const double minDensityFloor = 1e-6; // smaller floor — we will inspect weights
    const double maxWeight = 1e6; // essentially no cap here

    for (int E : energies) {
        printf("\n--- Energy %d GeV ---\n", E);
        std::string fname = fileOf[E];
        TFile f(fname.c_str());
        if (f.IsZombie()) { printf("  cannot open %s\n", fname.c_str()); continue; }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) { printf("  no Hits tree in %s\n", fname.c_str()); f.Close(); continue; }

        // branches
        Int_t detID, trackID, evID_branch;
        Double_t edep;
        t->SetBranchAddress("fDetectorID",&detID);
        t->SetBranchAddress("fEdep",&edep);
        t->SetBranchAddress("fTrackID",&trackID);
        t->SetBranchAddress("fEvent",&evID_branch);

        // first pass: per-event per-layer 7x7
        std::map<Int_t, std::vector<std::vector<std::vector<double>>>> perEventLayer;
        Long64_t nEntries = t->GetEntries();
        for (Long64_t i=0;i<nEntries;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            int L,cx,ry; decodeDet(detID,L,cx,ry);
            if (L<0 || L>=nLayers) continue;
            if (cx<0 || cx>=outSide || ry<0 || ry>=outSide) continue;
            if (perEventLayer.find(evID_branch) == perEventLayer.end())
                perEventLayer[evID_branch] = std::vector<std::vector<std::vector<double>>>(nLayers,
                    std::vector<std::vector<double>>(outSide, std::vector<double>(outSide,0.0)));
            perEventLayer[evID_branch][L][cx][ry] += edep;
        }
        f.Close();

        // compress to events: event id -> Event buffers (mapTot,mapPri,mapSec), compute centroid
        struct Ev {
            int peakL;
            std::vector<std::vector<double>> mapTot;
            std::vector<std::vector<double>> mapPri;
            std::vector<std::vector<double>> mapSec;
            double cx, cy;
            Ev():peakL(0),mapTot(eventSide,std::vector<double>(eventSide,0.0)),mapPri(eventSide,std::vector<double>(eventSide,0.0)),mapSec(eventSide,std::vector<double>(eventSide,0.0)),cx(0),cy(0){}
        };
        std::vector<std::pair<Int_t,Ev>> events; events.reserve(perEventLayer.size());

        for (auto &p : perEventLayer) {
            Int_t eid = p.first;
            auto &layerCell = p.second;
            std::vector<double> layerSum(nLayers,0.0);
            for (int L=0; L<nLayers; ++L) for (int x=0;x<outSide;++x) for (int y=0;y<outSide;++y) layerSum[L] += layerCell[L][x][y];
            int peakL = std::max_element(layerSum.begin(), layerSum.end()) - layerSum.begin();
            if (layerSum[peakL] <= 0) continue;
            Ev ev;
            ev.peakL = peakL;
            int offset = eventCenter - (outSide-1)/2;
            for (int x=0;x<outSide;++x) for (int y=0;y<outSide;++y) ev.mapTot[x+offset][y+offset] = layerCell[peakL][x][y];
            // centroid
            double sE=0,sX=0,sY=0;
            for (int x=0;x<eventSide;++x) for (int y=0;y<eventSide;++y) { double v=ev.mapTot[x][y]; sE+=v; sX+=v*x; sY+=v*y; }
            if (sE<=0) continue;
            ev.cx = sX/sE; ev.cy = sY/sE;
            events.emplace_back(eid, std::move(ev));
        }
        perEventLayer.clear();

        // build centroid density
        TH2D *hd = new TH2D("hd","centroid density", centroidBinsPerSide, -0.5, eventSide-0.5, centroidBinsPerSide, -0.5, eventSide-0.5);
        for (auto &pp : events) hd->Fill(pp.second.cx, pp.second.cy);
        // floor
        for (int ix=1; ix<=hd->GetNbinsX(); ++ix) for (int iy=1; iy<=hd->GetNbinsY(); ++iy) if (hd->GetBinContent(ix,iy) < minDensityFloor) hd->SetBinContent(ix,iy,minDensityFloor);

        // second pass: fill track-split maps for events (only for events vector). We reopen file and add pri/sec
        TFile f2(fname.c_str());
        if (f2.IsZombie()) { printf("  cannot re-open %s\n", fname.c_str()); delete hd; continue; }
        TTree *t2 = (TTree*)f2.Get("Hits");
        if (!t2) { printf("  no Hits tree second pass\n"); f2.Close(); delete hd; continue; }
        std::map<Int_t,int> evIndex;
        for (size_t i=0;i<events.size(); ++i) evIndex[events[i].first] = (int)i;
        Int_t det2, tr2, ev2; Double_t e2;
        t2->SetBranchAddress("fDetectorID",&det2);
        t2->SetBranchAddress("fEdep",&e2);
        t2->SetBranchAddress("fTrackID",&tr2);
        t2->SetBranchAddress("fEvent",&ev2);
        Long64_t n2 = t2->GetEntries();
        for (Long64_t i=0;i<n2;++i) {
            t2->GetEntry(i);
            if (e2 <= 0) continue;
            auto it = evIndex.find(ev2); if (it == evIndex.end()) continue;
            int idx = it->second;
            int L,cx,ry; decodeDet(det2,L,cx,ry);
            if (L != events[idx].second.peakL) continue;
            int offset = eventCenter - (outSide-1)/2;
            int ex = cx + offset, ey = ry + offset;
            if (ex<0||ex>=eventSide||ey<0||ey>=eventSide) continue;
            if (tr2 == 1) events[idx].second.mapPri[ex][ey] += e2;
            else events[idx].second.mapSec[ex][ey] += e2;
        }
        f2.Close();

        // final accumulators
        std::vector<std::vector<double>> accumTot(outSide, std::vector<double>(outSide,0.0));
        std::vector<std::vector<double>> accumPri(outSide, std::vector<double>(outSide,0.0));
        std::vector<std::vector<double>> accumSec(outSide, std::vector<double>(outSide,0.0));
        double sumW = 0.0;
        int nEvUsed = 0;

        for (size_t i=0;i<events.size(); ++i) {
            Ev &ed = events[i].second;
            // weight
            int bx = hd->GetXaxis()->FindBin(ed.cx);
            int by = hd->GetYaxis()->FindBin(ed.cy);
            double dens = hd->GetBinContent(bx,by);
            if (dens < minDensityFloor) dens = minDensityFloor;
            double w = makeFlat ? (1.0/dens) : 1.0;
            if (!TMath::Finite(w) || w<=0) w=1.0;
            if (w > maxWeight) w = maxWeight;
            // event total (for per-event normalize)
            double eTot=0;
            for (int x=0;x<eventSide;++x) for (int y=0;y<eventSide;++y) eTot += ed.mapTot[x][y];
            if (eTot<=0) continue;
            // shift to center
            int shiftX = 0, shiftY = 0;
            if (centerEvents) {
                shiftX = (int)TMath::Nint(eventCenter - ed.cx);
                shiftY = (int)TMath::Nint(eventCenter - ed.cy);
            }
            // crop outSide area
            for (int dx=-halfOut; dx<=halfOut; ++dx) for (int dy=-halfOut; dy<=halfOut; ++dy) {
                int sx = eventCenter + dx - shiftX;
                int sy = eventCenter + dy - shiftY;
                if (sx<0||sx>=eventSide||sy<0||sy>=eventSide) continue;
                int tx = dx + halfOut; int ty = dy + halfOut;
                double vT = perEventNormalize ? (ed.mapTot[sx][sy] / eTot) : ed.mapTot[sx][sy];
                double vP = perEventNormalize ? (ed.mapPri[sx][sy] / eTot) : ed.mapPri[sx][sy];
                double vS = perEventNormalize ? (ed.mapSec[sx][sy] / eTot) : ed.mapSec[sx][sy];
                accumTot[tx][ty] += w * vT;
                accumPri[tx][ty] += w * vP;
                accumSec[tx][ty] += w * vS;
            }
            sumW += w; nEvUsed++;
        }

        // divide by sumW to get average
        if (sumW <= 0) sumW = 1.0;
        // print arrays with high precision
        printf("Events used = %d   sumWeights = %.6f\n", nEvUsed, sumW);
        // helper to print and save CSV
        auto printAndSave = [&](const std::vector<std::vector<double>> &arr, const char* tag) {
            printf("\n--- %s (rows: y=0..6, cols: x=0..6) ---\n", tag);
            char outfn[256];
            snprintf(outfn, sizeof(outfn), "/tmp/map_E%d_%s.csv", E, tag);
            FILE *of = fopen(outfn,"w");
            if (!of) { printf("  cannot open %s for write\n", outfn); return; }
            // header
            fprintf(of,"col0,col1,col2,col3,col4,col5,col6\n");
            for (int y=0;y<outSide;++y) {
                for (int x=0;x<outSide;++x) {
                    double v = arr[x][y] / sumW;
                    printf("%12.8e ", v);
                    fprintf(of, "%.12e%s", v, (x==outSide-1) ? "" : "," );
                }
                printf("\n");
                fprintf(of, "\n");
            }
            fclose(of);
            printf("  saved /tmp/map_E%d_%s.csv\n", E, tag);
        };

        printAndSave(accumTot, "TOTAL");
        printAndSave(accumPri, "PRIMARY");
        printAndSave(accumSec, "SECONDARY");

        // simple radial profile print (avg at integer radius)
        int nRad = halfOut+1;
        std::vector<double> sumR(nRad,0.0);
        std::vector<int> cntR(nRad,0);
        for (int x=0;x<outSide;++x) for (int y=0;y<outSide;++y) {
            int dx = x - halfOut, dy = y - halfOut;
            double r = sqrt(dx*dx + dy*dy);
            int rb = (int)floor(r + 1e-9);
            if (rb<0) rb=0; if (rb>=nRad) rb=nRad-1;
            sumR[rb] += accumTot[x][y] / sumW;
            cntR[rb] += 1;
        }
        printf("\nRadial profile (radius, average):\n");
        for (int r=0;r<nRad;++r) {
            double val = (cntR[r]>0) ? sumR[r]/cntR[r] : 0.0;
            printf(" r=%d  %12.8e\n", r, val);
        }

        delete hd;
    } // energies loop

    printf("\nexportFinalMaps done. CSVs in /tmp/map_E<energy>_<TYPE>.csv\n");
}

