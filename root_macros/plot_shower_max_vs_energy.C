// plot_shower_max_vs_energy_converted_with_unitfix_and_logfit.C
// Same as previous macro but connects points and fits y = a*log(x) + b.
//
// Usage:
//  root -l plot_shower_max_vs_energy_converted_with_unitfix_and_logfit.C
//  plot_shower_max_vs_energy_converted_with_unitfix_and_logfit();  // default unitConv=1e-6

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

static std::map<std::pair<int,int>, double> readMuonCSV(const char* csv, double &overallMean) {
    std::ifstream ifs(csv);
    std::map<std::pair<int,int>, double> muMean;
    overallMean = 0;
    if (!ifs.is_open()) { printf("ERROR: cannot open mu CSV: %s\n", csv); return muMean; }
    std::string line; std::getline(ifs, line); // header
    double sum=0; int cnt=0;
    while (std::getline(ifs, line)) {
        if (line.size()<3) continue;
        std::stringstream ss(line);
        int col,row; double meanE; char comma;
        if (!(ss >> col)) continue; ss >> comma;
        if (!(ss >> row)) continue; ss >> comma;
        if (!(ss >> meanE)) continue;
        muMean[{col,row}] = meanE;
        sum += meanE; ++cnt;
    }
    ifs.close();
    if (cnt>0) overallMean = sum/double(cnt);
    return muMean;
}

static double readScale(const char* file) {
    std::ifstream ifs(file);
    if (!ifs.is_open()) return -1;
    std::string hdr,line; std::getline(ifs,hdr);
    if (!std::getline(ifs,line)) return -1;
    std::stringstream ss(line); double target,total,scale; char comma;
    if (!(ss>>target)) return -1; ss>>comma;
    if (!(ss>>total)) return -1; ss>>comma;
    if (!(ss>>scale)) return -1;
    return scale;
}

void plot_shower_max_vs_energy(
    const char* muCSV = "/home/rudradeb/sim/output/muon_mip_calibration.csv",
    const char* scaleFile = "/home/rudradeb/sim/output/electron_scale.txt",
    double unitConv = 1e-6 // multiply final peak energy by this before writing/plotting
)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    // edit these if needed
    std::vector<int> energies = {10,20,50,100};
    std::map<int,std::string> fileOf = {
        {10,  "/home/rudradeb/sim/muonoutput/output_thread-2.root"},
        {20,  "/home/rudradeb/sim/muonoutput/output_thread-3.root"},
        {50,  "/home/rudradeb/sim/muonoutput/output_thread-4.root"},
        {100, "/home/rudradeb/sim/muonoutput/output_thread-5.root"}
    };
    const int nLayers = 24;

    // read muon calibration
    double overallMean = 0;
    auto muMean = readMuonCSV(muCSV, overallMean);
    if (muMean.empty()) { printf("ERROR: muon CSV empty: %s\n", muCSV); return; }
    printf("Muon overall mean (units of fEdep) = %.6e\n", overallMean);

    // read global scale (GeV per MIP)
    double scale = readScale(scaleFile);
    if (scale <= 0) { printf("ERROR: bad scale file: %s\n", scaleFile); return; }
    printf("Global scale read: 1 MIP = %.12e GeV\n", scale);
    printf("Using unitConv = %.6e (multiply peak energy by this before output/plot)\n", unitConv);

    std::string outCSV = "/home/rudradeb/sim/output/shower_maxima_converted_fixed.csv";
    std::ofstream ofs(outCSV);
    if (!ofs.is_open()) { printf("ERROR: cannot write %s\n", outCSV.c_str()); return; }
    ofs << "EnergyGeV,nEvents,peakLayer,peakMeanEnergy_GeV_after_unitConv\n";

    std::vector<double> xE, yPeakEnergyScaled, yPeakLayer;

    for (int E : energies) {
        auto itf = fileOf.find(E);
        if (itf==fileOf.end()) { printf("No file for %d GeV\n", E); continue; }
        TFile f(itf->second.c_str());
        if (f.IsZombie()) { Warning("macro","Cannot open %s", itf->second.c_str()); continue; }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) { Warning("macro","No Hits tree in %s", itf->second.c_str()); f.Close(); continue; }

        Int_t detID=0; Double_t edep=0; Long64_t evID=-1; bool hasEvent=false;
        if (t->GetBranch("fDetectorID")) t->SetBranchAddress("fDetectorID",&detID);
        else if (t->GetBranch("detID")) t->SetBranchAddress("detID",&detID);
        if (t->GetBranch("fEdep")) t->SetBranchAddress("fEdep",&edep);
        else if (t->GetBranch("Edep")) t->SetBranchAddress("Edep",&edep);
        if (t->GetBranch("fEvent")) { t->SetBranchAddress("fEvent",&evID); hasEvent=true; }
        else if (t->GetBranch("event")) { t->SetBranchAddress("event",&evID); hasEvent=true; }
        else if (t->GetBranch("evt")) { t->SetBranchAddress("evt",&evID); hasEvent=true; }

        std::map<long long, std::vector<double>> evLayerGeV;
        Long64_t N = t->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            int L,cx,ry; decodeDet(detID,L,cx,ry);
            if (L<0 || L>=nLayers) continue;
            long long ev = hasEvent ? (long long)evID : (long long)i;
            if (!evLayerGeV.count(ev)) evLayerGeV[ev] = std::vector<double>(nLayers,0.0);

            double meanMip = overallMean;
            auto it = muMean.find({cx,ry});
            if (it != muMean.end() && it->second>0) meanMip = it->second;
            double nMips = (meanMip>0) ? (edep / meanMip) : 0.0;
            double eGeV = nMips * scale;
            evLayerGeV[ev][L] += eGeV;
        }
        f.Close();

        int nEvents = evLayerGeV.size();
        if (nEvents==0) { Warning("macro","No events for %d GeV", E); continue; }

        std::vector<double> meanLayer(nLayers,0.0);
        for (auto &p: evLayerGeV) for (int L=0; L<nLayers; ++L) meanLayer[L] += p.second[L];
        for (int L=0; L<nLayers; ++L) meanLayer[L] /= double(nEvents);

        int peakL = std::distance(meanLayer.begin(), std::max_element(meanLayer.begin(), meanLayer.end()));
        double peakMeanGeV = meanLayer[peakL];

        // apply user conversion factor (unitConv)
        double peakMeanGeV_scaled = peakMeanGeV * unitConv;

        printf("E=%d GeV: nEvents=%d peakL=%d raw_peak=%.12e  scaled_peak=%.12e\n",
               E, nEvents, peakL, peakMeanGeV, peakMeanGeV_scaled);

        ofs << E << "," << nEvents << "," << peakL << "," << std::setprecision(12) << peakMeanGeV_scaled << "\n";

        xE.push_back(double(E));
        yPeakEnergyScaled.push_back(peakMeanGeV_scaled);
        yPeakLayer.push_back(double(peakL));
    }

    ofs.close();
    printf("Wrote %s\n", outCSV.c_str());

    if (xE.empty()) { printf("No points to plot.\n"); return; }

    // plot
    TCanvas *c = new TCanvas("c_peakE","Peak energy (scaled) vs beam",900,600);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
    gPad->SetLogx(1);

    double xmin = *std::min_element(xE.begin(),xE.end());
    double xmax = *std::max_element(xE.begin(),xE.end());
    double ymin = *std::min_element(yPeakEnergyScaled.begin(), yPeakEnergyScaled.end());
    double ymax = *std::max_element(yPeakEnergyScaled.begin(), yPeakEnergyScaled.end());

    TH1D *frame = new TH1D("frame","Peak mean energy (scaled);Beam energy (GeV);Peak mean energy (scaled units)",
                           10, xmin*0.8, xmax*1.4);
    frame->SetMinimum(std::max(0.0, ymin*0.8));
    frame->SetMaximum(ymax*1.4);
    frame->Draw();

    // Graph with connected points
    TGraph *g = new TGraph(xE.size());
    for (size_t i=0;i<xE.size();++i) g->SetPoint(i, xE[i], yPeakEnergyScaled[i]);
    g->SetMarkerStyle(21); g->SetMarkerSize(1.4); g->SetLineWidth(2);
    g->SetLineColor(kBlue);
    g->Draw("LP SAME"); // points connected by a line

    // logarithmic fit: y = a * log(x) + b
    TF1 *flog = new TF1("flog","[0]*log(x)+[1]", xmin*0.9, xmax*1.1);
    flog->SetParNames("a","b");
    g->Fit(flog,"RQ"); // R: use fit range, Q: quiet

    // draw fitted curve (dense)
    flog->SetLineColor(kRed);
    flog->SetLineWidth(2);
    flog->Draw("LSAME");

    // annotate fit parameters and chi2/ndf
    double a = flog->GetParameter(0);
    double b = flog->GetParameter(1);
    double chi2 = flog->GetChisquare();
    double ndf  = flog->GetNDF();
    TLatex t; t.SetNDC(); t.SetTextSize(0.035);
    t.DrawLatex(0.15,0.85,Form("fit: y = a ln(x) + b"));
    t.DrawLatex(0.15,0.80,Form("a = %.5g", a));
    t.DrawLatex(0.15,0.75,Form("b = %.5g", b));
    t.DrawLatex(0.15,0.70,Form("#chi^{2}/ndf = %.2f / %.0f", chi2, ndf));

    // legend
    TLegend *leg = new TLegend(0.62,0.65,0.92,0.90);
    leg->SetFillColor(0); leg->SetBorderSize(0);
    leg->AddEntry(g, "data (connected)", "lp");
    leg->AddEntry(flog, "log fit", "l");
    leg->Draw();

    // label each point with peak layer and value
    for (size_t i=0;i<xE.size();++i) {
        TLatex label; label.SetTextSize(0.03);
        label.DrawLatex(xE[i]*1.02, yPeakEnergyScaled[i]*1.03,
                        Form("L=%d, %.3g", int(yPeakLayer[i]), yPeakEnergyScaled[i]));
    }

    c->SaveAs("/home/rudradeb/sim/output/peakEnergy_vs_beam_fixed_logfit.png");
    printf("Saved /home/rudradeb/sim/output/peakEnergy_vs_beam_fixed_logfit.png\n");

    // cleanup
    delete frame;
    delete g;
    delete flog;
    delete leg;
    delete c;
}

