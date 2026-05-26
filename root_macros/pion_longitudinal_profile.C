#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

void pion_longitudinal_profile() {
    gStyle->SetOptStat(0);

    const char* fileName = "/home/rudradeb/sim/output/output_thread-1.root";
    const int nLayers = 24;

    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    TTree* tree = (TTree*)f->Get("Hits");
    if (!tree) {
        std::cerr << "Tree Hits not found\n";
        return;
    }

    Int_t detID;
    Double_t edep;
    tree->SetBranchAddress("fDetectorID", &detID);
    tree->SetBranchAddress("fEdep", &edep);

    std::vector<double> layerSum(nLayers, 0.0);

    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);

        int layer = (detID / 10000) - 100;   // active detIDs have (layer+100)*10000
        if (layer >= 0 && layer < nLayers)
            layerSum[layer] += edep;
    }

    TCanvas* c = new TCanvas("c", "Pion Longitudinal Shower Profile", 800, 600);
    TGraph* g = new TGraph(nLayers);

    for (int l = 0; l < nLayers; ++l)
        g->SetPoint(l, l + 1, layerSum[l]);

    g->SetTitle("Pion Longitudinal Shower Profile;Layer;Deposited Energy (MeV)");
    g->SetMarkerStyle(20);
    g->SetLineWidth(2);
    g->Draw("APL");

    c->SaveAs("pion_longitudinal_profile.png");
}

