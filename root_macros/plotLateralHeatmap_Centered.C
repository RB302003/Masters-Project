// plot_lateral_heatmaps_fixed.C
// Center per-event showers using entrance-layer(s), then sum peak-layer energy.
// Usage: plotLateralHeatmap_CenteredFixed(100, kFALSE, 2)
// energy values supported: 10,20,50,100 (change file paths below to match your files)

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMarker.h"

#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
using std::cout;
using std::endl;

void plotLateralHeatmap_Centered(int energy = 100, Bool_t normalize = kFALSE,
                                 int useFirstNLayersForImpact = 1 /* try 1 or 2 */) {
  // aesthetic
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  // map energy -> file (update these paths if needed)
  std::map<int,std::string> fileOf = {
    {10,  "/home/rudradeb/sim/finaloutput/output_thread-1.root"},
    {20,  "/home/rudradeb/sim/finaloutput/output_thread-2.root"},
    {50,  "/home/rudradeb/sim/finaloutput/output_thread-3.root"},
    {100, "/home/rudradeb/sim/finaloutput/output_thread-4.root"}
  };
  if (!fileOf.count(energy)) {
    printf("Unknown energy %d\n", energy);
    return;
  }
  const std::string fname = fileOf[energy];

  // detector grid parameters (7x7)
  const int nSide = 7;
  const int colMin = 0, colMax = nSide-1;
  const int rowMin = 0, rowMax = nSide-1;
  // chosen center grid indices (0..6). For 7x7 center is (3,3)
  const int cCol = 3, cRow = 3;

  // decode packed detID (layer*10000 + col*100 + row)
  auto decode = [](int detID, int &layer, int &col, int &row){
    layer = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int local = detID % 10000;
    col = local / 100;
    row = local % 100;
  };

  // --- Pass 0: open file and determine which layer indices exist, and find peak layer ---
  std::set<int> layersSeen;
  std::map<int,double> layerSum; // layer -> total energy
  {
    TFile f(fname.c_str());
    if (f.IsZombie()) { printf("Cannot open %s\n", fname.c_str()); return; }
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) { printf("No Hits tree in %s\n", fname.c_str()); f.Close(); return; }

    Int_t detID; Double_t edep;
    t->SetBranchAddress("fDetectorID",&detID);
    t->SetBranchAddress("fEdep",&edep);

    Long64_t n = t->GetEntries();
    for (Long64_t i=0;i<n;++i) {
      t->GetEntry(i);
      if (edep <= 0) continue;
      int L,cx,ry; decode(detID,L,cx,ry);
      layersSeen.insert(L);
      layerSum[L] += edep;
    }
    f.Close();
  }

  if (layersSeen.empty()) { printf("No layers found in file %s\n", fname.c_str()); return; }

  // sort layers and pick lowest N as "entrance" layers
  std::vector<int> layersVec(layersSeen.begin(), layersSeen.end());
  std::sort(layersVec.begin(), layersVec.end());
  int availableLayers = (int)layersVec.size();
  printf("Found %d distinct layers (lowest..highest):", availableLayers);
  for (int L : layersVec) printf(" %d", L);
  printf("\n");

  // decide entrance layers to use: lowest useFirstNLayersForImpact layers
  int useN = std::max(1, useFirstNLayersForImpact);
  useN = std::min(useN, (int)layersVec.size());
  std::vector<int> entranceLayers(layersVec.begin(), layersVec.begin()+useN);
  printf("Using entrance layers (for impact detection):");
  for (int L : entranceLayers) printf(" %d", L);
  printf("\n");

  // find peak layer by total energy
  int peakLayer = -1;
  {
    double best = -1;
    for (auto &kv : layerSum) {
      if (kv.second > best) { best = kv.second; peakLayer = kv.first; }
    }
    printf("Peak layer chosen = %d (total E = %g)\n", peakLayer, layerSum[peakLayer]);
  }

  // --- Pass 1: per-event entrance sums to infer impact cell ---
  // builds map: eventID -> (col,row) of impact cell
  std::map<int,std::pair<int,int>> impactCellOfEvent;
  {
    TFile f(fname.c_str());
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) { f.Close(); return; }

    Int_t detID, eventID; Double_t edep;
    t->SetBranchAddress("fDetectorID",&detID);
    t->SetBranchAddress("fEdep",&edep);
    t->SetBranchAddress("fEvent",&eventID);

    // temporary per-event 7x7 grids for entrance layers
    std::map<int, std::vector<std::vector<double>>> entrSum; // event -> grid

    Long64_t n = t->GetEntries();
    for (Long64_t i=0;i<n;++i) {
      t->GetEntry(i);
      if (edep <= 0) continue;
      int L,cx,ry; decode(detID,L,cx,ry);
      // only count hits in chosen entrance layers
      if (std::find(entranceLayers.begin(), entranceLayers.end(), L) == entranceLayers.end()) continue;
      if (cx < colMin || cx > colMax || ry < rowMin || ry > rowMax) continue;
      auto &grid = entrSum[eventID];
      if (grid.empty()) grid.assign(nSide, std::vector<double>(nSide, 0.0));
      grid[cx][ry] += edep;
    }

    // pick best cell per event
    for (auto &kv : entrSum) {
      int ev = kv.first;
      auto &grid = kv.second;
      int bestC = cCol, bestR = cRow;
      double bestE = -1;
      for (int cx = colMin; cx <= colMax; ++cx) {
        for (int ry = rowMin; ry <= rowMax; ++ry) {
          if (grid[cx][ry] > bestE) { bestE = grid[cx][ry]; bestC = cx; bestR = ry; }
        }
      }
      impactCellOfEvent[ev] = {bestC, bestR};
    }
    printf("Found impact cells from entrance layers for %zu events\n", impactCellOfEvent.size());
    f.Close();
  }

  // --- Pass 2: accumulate peak-layer hits, recentered per event ---
  std::vector<std::vector<double>> acc(nSide, std::vector<double>(nSide, 0.0));
  int usedEvents = 0;
  long long totalPeakHits = 0;

  {
    TFile f(fname.c_str());
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) { f.Close(); return; }

    Int_t detID, eventID; Double_t edep;
    t->SetBranchAddress("fDetectorID",&detID);
    t->SetBranchAddress("fEdep",&edep);
    t->SetBranchAddress("fEvent",&eventID);

    Long64_t n = t->GetEntries();
    // We'll also build per-event peak-layer energy to optionally fallback to peak-cell per-event
    std::map<int, std::vector<std::vector<double>>> peakGridPerEvent;

    for (Long64_t i=0;i<n;++i) {
      t->GetEntry(i);
      if (edep <= 0) continue;
      int L,cx,ry; decode(detID,L,cx,ry);
      if (L != peakLayer) continue;
      if (cx < colMin || cx > colMax || ry < rowMin || ry > rowMax) continue;

      peakGridPerEvent[eventID]; // ensure entry exists
      peakGridPerEvent[eventID].empty();
      auto &grid = peakGridPerEvent[eventID];
      if (grid.empty()) grid.assign(nSide, std::vector<double>(nSide, 0.0));
      grid[cx][ry] += edep;
      totalPeakHits++;
    }

    // Now for each event with peak-layer data, find impact cell (entrance) if present; otherwise fallback to peak-layer max cell
    for (auto &kv : peakGridPerEvent) {
      int ev = kv.first;
      auto &pgrid = kv.second;

      // determine this event impact cell
      int impC = cCol, impR = cRow;
      auto itImp = impactCellOfEvent.find(ev);
      if (itImp != impactCellOfEvent.end()) {
        impC = itImp->second.first;
        impR = itImp->second.second;
      } else {
        // fallback: use this event's strongest peak-layer cell as impact (keeps event instead of discarding)
        double bestE = -1;
        for (int cx = colMin; cx <= colMax; ++cx)
          for (int ry = rowMin; ry <= rowMax; ++ry)
            if (pgrid[cx][ry] > bestE) { bestE = pgrid[cx][ry]; impC = cx; impR = ry; }
      }

      // now shift and add pgrid into acc
      int sx = cCol - impC;
      int sy = cRow - impR;

      bool anyAdded = false;
      for (int cx = colMin; cx <= colMax; ++cx) {
        for (int ry = rowMin; ry <= rowMax; ++ry) {
          int cx2 = cx + sx;
          int ry2 = ry + sy;
          if (cx2 < colMin || cx2 > colMax || ry2 < rowMin || ry2 > rowMax) continue;
          double val = pgrid[cx][ry];
          if (val <= 0) continue;
          acc[cx2][ry2] += val;
          anyAdded = true;
        }
      }
      if (anyAdded) ++usedEvents;
    }

    f.Close();
  } // end pass 2

  printf("Total peak-layer hit entries processed: %lld\n", totalPeakHits);
  printf("Number of events contributing to final centered map: %d\n", usedEvents);

  // --- Draw heatmap ---
  TCanvas *c = new TCanvas("c_centered","Centered lateral heatmap",1000,850);

  // convert to TH2D. Bins centered at 0..6 (columns/rows)
  TH2D *h = new TH2D("h_center","Centered Lateral Cross-Section;Column (cell index);Row (cell index)",
                     nSide, -0.5, nSide-0.5, nSide, -0.5, nSide-0.5);

  double maxVal = 0;
  for (int x=0;x<nSide;++x) for (int y=0;y<nSide;++y) maxVal = std::max(maxVal, acc[x][y]);
  if (normalize && maxVal <= 0) maxVal = 1.0;

  // Fill bins. Use bin numbers directly to avoid ambiguity.
  for (int x=0;x<nSide;++x) {
    for (int y=0;y<nSide;++y) {
      double v = normalize ? (acc[x][y] / maxVal) : acc[x][y];
      int binx = h->GetXaxis()->FindBin(x); // FindBin on coordinate x (0..6) returns bin number
      int biny = h->GetYaxis()->FindBin(y);
      h->SetBinContent(binx, biny, v);
    }
  }

  h->SetContour(100);
  h->Draw("COLZ");

  TLatex lat; lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.14,0.94, Form("Centered lateral cross-section @ peak layer %d  (%d GeV)", peakLayer, energy));
  lat.SetTextSize(0.03);
  lat.DrawLatex(0.14,0.90, normalize ? "Normalized to max bin (centered)" : "Absolute energy (sum, centered)");

  // mark center cell
  TMarker m(cCol, cRow, 29);
  m.SetMarkerColor(kWhite);
  m.SetMarkerSize(1.6);
  m.Draw();

  c->Update();
}

