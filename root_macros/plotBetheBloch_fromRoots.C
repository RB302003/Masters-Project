// plotBetheBloch_fromRoots_fixed.C
// Safer variant with extra checks and no unsafe format strings.
// Usage:
//   root -l plotBetheBloch_fromRoots_fixed.C
//   then: plotBetheBloch_fromRoots_fixed();

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

static inline void decodeDet(int detID, int &L, int &col, int &row){
    L  = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

void plotBetheBloch_fromRoots(bool useLayerRange = false,
                              int layerSingle = 4,
                              int layerFrom = 2, int layerTo = 6,
                              double energyScale = 1e-3)
{
    // === EDIT: put your actual files here ===
    std::map<double,std::string> files = {
        {10.0, "/home/rudradeb/sim/muonoutput/output_thread-1.root"},
        {20.0, "/home/rudradeb/sim/muonoutput/output_thread-2.root"},
        {50.0, "/home/rudradeb/sim/muonoutput/output_thread-3.root"},
        {100.0,"/home/rudradeb/sim/muonoutput/output_thread-4.root"}
    };
    // ========================================

    std::vector<double> energies;
    std::vector<double> mpvs, mpvErrs;

    for (auto &kv : files) {
        double Ebeam = kv.first;
        std::string fname = kv.second;
        TFile fin(fname.c_str());
        if (!fin.IsOpen() || fin.IsZombie()) {
            Warning("plotBetheBloch_fromRoots_fixed","Cannot open %s -- skipping", fname.c_str());
            continue;
        }
        TTree *tree = (TTree*) fin.Get("Hits");
        if (!tree) {
            Warning("plotBetheBloch_fromRoots_fixed","No Hits TTree in %s -- skipping", fname.c_str());
            fin.Close(); continue;
        }

        Int_t detID=0;
        Double_t edep=0;
        tree->SetBranchAddress("fDetectorID",&detID);
        tree->SetBranchAddress("fEdep",&edep);

        // histogram parameters (units: GeV after energyScale)
        double hmin = 0.0, hmax = 0.5;
        TH1D *h = new TH1D("h_tmp","Edep per layer",200,hmin,hmax);

        Long64_t N = tree->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            tree->GetEntry(i);
            if (edep <= 0) continue;
            int L,col,row; decodeDet(detID,L,col,row);
            if (useLayerRange) {
                if (L < layerFrom || L > layerTo) continue;
            } else {
                if (L != layerSingle) continue;
            }
            h->Fill(edep * energyScale);
        }

        Long64_t nEntries = h->GetEntries();
        if (nEntries < 10) {
            Warning("plotBetheBloch_fromRoots_fixed",
                    "Low statistics for %.1f GeV (%lld entries) -- skipping", Ebeam, (long long)nEntries);
            delete h;
            fin.Close();
            continue;
        }

        // Fit Landau with quiet options, and check fit return status
        TF1 *fl = new TF1("fl_tmp","landau", hmin, hmax);
        fl->SetParameters(h->GetMaximum(), h->GetBinCenter(h->GetMaximumBin()), 0.02);
        int fitStatus = h->Fit(fl, "RQ"); // R=use fit range, Q=quiet
        double mpv = fl->GetParameter(1);
        double mpvErr = fl->GetParError(1);

        if (fitStatus != 0 || !std::isfinite(mpv) || mpv <= hmin || mpv >= hmax) {
            // fallback to histogram mode (peak bin) if fit fails
            mpv = h->GetBinCenter(h->GetMaximumBin());
            mpvErr = h->GetRMS() / sqrt(std::max(1.0, double(nEntries)));
            Warning("plotBetheBloch_fromRoots_fixed","Fit failed or unreliable for %.1f GeV: using mode (%.4g +/- %.4g)",
                    Ebeam, mpv, mpvErr);
        }

        energies.push_back(Ebeam);
        mpvs.push_back(mpv);
        mpvErrs.push_back(mpvErr);

        delete fl;
        delete h;
        fin.Close();
    }

    if (energies.empty()) {
        Error("plotBetheBloch_fromRoots_fixed","No valid input samples processed. Check file list and permissions.");
        return;
    }

    // after you have filled energies[], mpvs[], mpvErrs[] and created g (TGraphErrors)
int n = energies.size();
TGraphErrors *g = new TGraphErrors(n);
for (int i=0;i<n;++i) {
    g->SetPoint(i, energies[i], mpvs[i]);
    g->SetPointError(i, 0.0, mpvErrs[i]);
}

// Styling
g->SetMarkerStyle(21);
g->SetMarkerSize(1.0);
g->SetLineWidth(2);
g->SetLineColor(kBlack);
g->SetMarkerColor(kBlack);
g->SetTitle("MPV of per-layer E_{dep} vs beam energy;Beam energy (GeV);MPV E_{dep} per layer (GeV)");

// Canvas: log x
TCanvas *c = new TCanvas("cBB","Bethe-Bloch style",900,650);
c->SetLogx();                 // log x-axis
gPad->SetLeftMargin(0.12);
gPad->SetBottomMargin(0.12);

// Draw points and straight connecting line
g->Draw("ALP");               // A axes L line P points

// Draw a smooth spline through the points for visual guidance
// Use TSpline3 which will interpolate; the spline respects the graph x-values even on log x
TSpline3 *s = new TSpline3("spline", g);
s->SetLineWidth(2);
s->SetLineColor(kBlue);
s->SetLineStyle(1);
s->Draw("same");              // overlay

// Optionally draw a faint connecting polyline with smoothing (for emphasis)
TGraph *gLine = new TGraph(n);
for (int i=0;i<n;++i) gLine->SetPoint(i, energies[i], mpvs[i]);
gLine->SetLineColor(kGray+2);
gLine->SetLineStyle(7);       // dashed
gLine->SetLineWidth(1);
gLine->Draw("Lsame");

// Add legend
TLegend *leg = new TLegend(0.58,0.75,0.88,0.88);
leg->SetBorderSize(0);
leg->AddEntry(g, "MPV (data)", "lp");
leg->AddEntry(s, "Spline guide", "l");
leg->Draw();

// redraw to update axis after log scale
c->Update();
c->SaveAs("MPV_vs_energy_smoothed.png");

}

