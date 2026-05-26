// generate_MPV_plot_fixed.C
// Usage:
//   root -l generate_MPV_plot_fixed.C
//   then at the ROOT prompt: generate_MPV_plot()

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
#include <algorithm>

static inline void decodeDet(int detID, int &L, int &col, int &row){
    L  = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

void generate_MPV_plot(bool useLayerRange = true,
                       int layerSingle = 4,
                       int layerFrom = 2, int layerTo = 6,
                       double energyScale = 1e-3, // stored units -> GeV (adjust if needed)
                       int nbins = 200,
                       double hmin = 0.0, double hmax = 0.5,
                       const char *outRootDebug = "mpv_debug.root",
                       const char *outPNG = "MPV_vs_energy.png")
{
    // --- EDIT THIS MAP: energy(GeV) -> file path ---
    std::map<double,std::string> files = {
        {10.0, "/home/rudradeb/sim/muonoutput/output_thread-4.root"},
        {20.0, "/home/rudradeb/sim/muonoutput/output_thread-5.root"},
        {50.0, "/home/rudradeb/sim/muonoutput/output_thread-6.root"},
        {100.0,"/home/rudradeb/sim/muonoutput/output_thread-7.root"}
    };
    // ------------------------------------------------

    std::vector<double> energies;
    std::vector<double> mpvs;
    std::vector<double> mpvErrs;
    TFile fout(outRootDebug, "RECREATE");

    for (auto &kv : files) {
        double Ebeam = kv.first;
        std::string fname = kv.second;
        TFile fin(fname.c_str());
        if (!fin.IsOpen() || fin.IsZombie()) {
            Warning("generate_MPV_plot","Cannot open %s -- skipping", fname.c_str());
            continue;
        }
        TTree *tree = dynamic_cast<TTree*>(fin.Get("Hits"));
        if (!tree) {
            Warning("generate_MPV_plot","No Hits TTree in %s -- skipping", fname.c_str());
            fin.Close();
            continue;
        }

        Int_t detID = 0;
        Double_t edep = 0;
        tree->SetBranchAddress("fDetectorID",&detID);
        tree->SetBranchAddress("fEdep",&edep);

        // choose histogram name per energy to avoid collisions
        std::string hname = Form("hEdep_%.0fGeV", Ebeam);
        TH1D *h = new TH1D(hname.c_str(), Form("Edep per hit (%.0f GeV);E_{dep} (GeV);Counts", Ebeam),
                           nbins, hmin, hmax);

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

        Long64_t entries = h->GetEntries();
        if (entries < 20) {
            Warning("generate_MPV_plot","Low statistics for %.0f GeV (%lld entries) -- skipping", Ebeam, (long long)entries);
            delete h;
            fin.Close();
            continue;
        }

        // Fit with Landau; fall back to mode if fit fails
        TF1 *fl = new TF1(Form("fl_%.0f", Ebeam),"landau", hmin, hmax);
        fl->SetParameters(h->GetMaximum(), h->GetBinCenter(h->GetMaximumBin()), (hmax-hmin)/50.0);
        int fitStatus = h->Fit(fl, "RQ"); // R=range, Q=quiet
        double mpv = fl->GetParameter(1);
        double mpvErr = fl->GetParError(1);

        if (fitStatus != 0 || !std::isfinite(mpv) || mpv <= hmin || mpv >= hmax) {
            // fallback: histogram mode + estimated error
            mpv = h->GetBinCenter(h->GetMaximumBin());
            mpvErr = h->GetRMS() / sqrt(std::max(1.0, double(entries)));
            Warning("generate_MPV_plot","Fit failed/unreliable for %.0f GeV: using mode (%.4g +/- %.4g)", Ebeam, mpv, mpvErr);
        }

        // Save histogram and fit to debug ROOT file
        h->Write();
        fl->Write(Form("fit_%.0f", Ebeam));

        energies.push_back(Ebeam);
        mpvs.push_back(mpv);
        mpvErrs.push_back(mpvErr);

        delete fl;
        delete h;
        fin.Close();
    }

    // done with input loop
    if (energies.empty()) {
        Error("generate_MPV_plot","No valid runs processed. Edit file map or check files.");
        fout.Close();
        return;
    }

    // build graph and draw
    int n = energies.size();
    TGraphErrors *g = new TGraphErrors(n);
    for (int i=0;i<n;++i) {
        g->SetPoint(i, energies[i], mpvs[i]);
        g->SetPointError(i, 0.0, mpvErrs[i]);
    }
    g->SetMarkerStyle(21); g->SetMarkerSize(1.2);
    g->SetLineWidth(2);
    g->SetTitle("MPV of per-layer E_{dep} vs beam energy;Beam energy (GeV);MPV E_{dep} per layer (GeV)");

    // Canvas (log x)
    TCanvas *c = new TCanvas("cMPV","MPV vs Energy",900,650);
    c->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    g->Draw("AP");

    // connect points by building a simple TGraph from the vectors
    TGraph *gline = new TGraph(n);
    for (int i=0;i<n;++i) gline->SetPoint(i, energies[i], mpvs[i]);
    gline->SetLineStyle(1); gline->SetLineWidth(1); gline->SetLineColor(kGray+2);
    gline->Draw("Lsame");

    // compute energy range for log-fit safely
    double emin = *std::min_element(energies.begin(), energies.end());
    double emax = *std::max_element(energies.begin(), energies.end());
    if (emin <= 0) emin = 1e-6;
    double fit_lo = 0.9 * emin;
    double fit_hi = 1.1 * emax;

    // optionally fit a log model: y = a + b*ln(E)
    TF1 *f_log = new TF1("f_log","[0] + [1]*log(x)", fit_lo, fit_hi);
    f_log->SetParameters(mpvs.front(), 1e-3);
    int stat = g->Fit(f_log, "R"); // fit in range
    if (stat == 0) {
        f_log->SetLineColor(kRed);
        f_log->SetLineWidth(2);
        f_log->Draw("same");
    } else {
        Warning("generate_MPV_plot","Log fit failed or unreliable (status=%d).", stat);
    }

    // Legend
    TLegend *leg = new TLegend(0.58,0.75,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(g,"MPV (data)","lp");
    if (stat==0) leg->AddEntry(f_log,"a + b ln(E)","l");
    leg->Draw();

    c->SaveAs(outPNG);

    // Write graph and fit to debug file
    g->Write("gMPV");
    if (f_log) f_log->Write("f_log");

    // cleanup
    delete gline;
    delete g;
    delete f_log;
    delete c;

    fout.Close();

    std::cout << "Done. Saved " << outPNG << " and debug file " << outRootDebug << std::endl;
}

