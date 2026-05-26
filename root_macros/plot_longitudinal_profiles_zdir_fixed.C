#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <iostream>
#include <vector>

void plot_longitudinal_profiles_zdir(const char* filename = "finaloutput/output_thread-3.root") {
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file!\n";
        return;
    }

    TTree *tree = (TTree*)file->Get("Hits");
    if (!tree) {
        std::cerr << "Hits tree not found!\n";
        return;
    }

    double edep;
    int detid;
    tree->SetBranchAddress("fEdep", &edep);
    tree->SetBranchAddress("fDetectorID", &detid);

    const int nLayers = 24;
    const int nRows = 7;  
    const int nCols = 7;

    std::vector<std::vector<double>> sliceProfile(nRows, std::vector<double>(nLayers, 0));
    std::vector<std::vector<double>> colProfile(nCols, std::vector<double>(nLayers, 0));

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (edep <= 0) continue;

        int layer = (detid / 10000) - 100;  // active detIDs have (layer+100)*10000
        int col   = (detid % 10000) / 100;   // 0 to 6
        int row   = (detid % 10000) % 100;   // 0 to 6

        if (layer >= 0 && layer < nLayers && row >= 0 && row < nRows && col >= 0 && col < nCols) {
            sliceProfile[row][layer] += edep;
            colProfile[col][layer] += edep;
        }
    }

    // --- Plot: Longitudinal by Row (slice) ---
    TCanvas *c1 = new TCanvas("c1", "Longitudinal Profiles by Slice (Row)", 800, 600);
    TMultiGraph *mgSlice = new TMultiGraph();
    TLegend *leg1 = new TLegend(0.7, 0.6, 0.9, 0.9);

    for (int r = 0; r < nRows; ++r) {
        TGraph *g = new TGraph();
        for (int l = 0; l < nLayers; ++l)
            g->SetPoint(l, l + 1, sliceProfile[r][l]);

        g->SetMarkerStyle(20);
        g->SetMarkerSize(1);
        g->SetLineWidth(2);
        g->SetLineColor(r + 1);
        g->SetMarkerColor(r + 1);

        mgSlice->Add(g);
        leg1->AddEntry(g, Form("Row %d", r), "lp");
    }

    mgSlice->SetTitle("Longitudinal Profiles by Row;Layer (Z);Total Edep (MeV)");
    mgSlice->Draw("APL");
    leg1->Draw();
    c1->SaveAs("longitudinal_by_row_3D.pdf");

    // --- Plot: Longitudinal by Column ---
    TCanvas *c2 = new TCanvas("c2", "Longitudinal Profiles by Column", 800, 600);
    TMultiGraph *mgCol = new TMultiGraph();
    TLegend *leg2 = new TLegend(0.7, 0.6, 0.9, 0.9);

    for (int c = 0; c < nCols; ++c) {
        TGraph *g = new TGraph();
        for (int l = 0; l < nLayers; ++l)
            g->SetPoint(l, l + 1, colProfile[c][l]);

        g->SetMarkerStyle(21);
        g->SetMarkerSize(1);
        g->SetLineWidth(2);
        g->SetLineColor(c + 1);
        g->SetMarkerColor(c + 1);

        mgCol->Add(g);
        leg2->AddEntry(g, Form("Col %d", c), "lp");
    }

    mgCol->SetTitle("Longitudinal Profiles by Column;Layer (Z);Total Edep (MeV)");
    mgCol->Draw("APL");
    leg2->Draw();
    c2->SaveAs("longitudinal_by_column_3D.pdf");

    file->Close();
}

