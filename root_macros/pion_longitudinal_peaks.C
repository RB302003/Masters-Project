#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <vector>
#include <iostream>

void pion_longitudinal_peaks() {
    gStyle->SetOptStat(0);

    const char* fileName = "/home/rudradeb/sim/output/output_thread-1.root";
    const int nLayers = 24;

    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) return;

    TTree* tree = (TTree*)f->Get("Hits");
    if (!tree) return;

    Int_t detID;
    Double_t edep;
    tree->SetBranchAddress("fDetectorID", &detID);
    tree->SetBranchAddress("fEdep", &edep);

    std::vector<double> peakE(nLayers, 0.0);

    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);

        int layer = (detID / 10000) - 100;
        if (layer >= 0 && layer < nLayers)
            if (edep > peakE[layer]) peakE[layer] = edep;
    }

    TCanvas* c = new TCanvas("c", "Peak Energy per Layer (Pions)", 800, 600);
    TGraph* g = new TGraph(nLayers);

    for (int l = 0; l < nLayers; ++l)
        g->SetPoint(l, l + 1, peakE[l]);

    g->SetTitle("Pion Shower: Peak Energy per Layer;Layer;Peak Energy (MeV)");
    g->SetMarkerStyle(21);
    g->SetLineWidth(2);
    g->Draw("APL");

    c->SaveAs("pion_longitudinal_peaks.png");
}

