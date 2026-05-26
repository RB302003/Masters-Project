// peakEnergy_converted_vs_beam.C
// Convert peak-layer mean energy into GeV using muon per-cell calibration + global scale,
// then plot peak (mean energy in peak layer) vs beam energy.
//
// Usage:
//   root -l peakEnergy_converted_vs_beam.C
//   peakEnergy_converted_vs_beam();

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLatex.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

// read muon CSV created by muon_calibrate: header then lines "col,row,meanE,calibFactor,...".
// Returns map[(col,row)] = meanMIP (mean energy deposit for a MIP in that cell).
static std::map<std::pair<int,int>, double> readMuonCSV(const char* muCSV, double &overallMean) {
    std::map<std::pair<int,int>, double> muMean;
    overallMean = 0.0;
    std::ifstream ifs(muCSV);
    if (!ifs.is_open()) {
        printf("WARNING: cannot open muon CSV '%s'\n", muCSV);
        return muMean;
    }
    std::string line;
    // skip header
    if (!std::getline(ifs, line)) { ifs.close(); return muMean; }
    double sum = 0.0;
    int cnt = 0;
    while (std::getline(ifs, line)) {
        if (line.size() == 0) continue;
        std::stringstream ss(line);
        int col, row;
        double meanE = 0.0;
        // try parse "col,row,meanE,..." robustly
        char comma;
        if (!(ss >> col)) continue;
        ss >> comma;
        if (!(ss >> row)) continue;
        ss >> comma;
        if (!(ss >> meanE)) meanE = 0.0;
        muMean[{col,row}] = meanE;
        if (meanE > 0) { sum += meanE; ++cnt; }
    }
    ifs.close();
    if (cnt>0) overallMean = sum / double(cnt);
    else overallMean = 0.0;
    return muMean;
}

// read simple scale file: header then line "targetEnergy,total_MIPs,scale_GeV_per_MIP"
static double readScale(const char* scaleFile) {
    std::ifstream ifs(scaleFile);
    if (!ifs.is_open()) { printf("WARNING: cannot open scale file '%s'\n", scaleFile); return -1.0; }
    std::string line;
    if (!std::getline(ifs, line)) { ifs.close(); return -1.0; } // header
    if (!std::getline(ifs, line)) { ifs.close(); return -1.0; }
    std::stringstream ss(line);
    double target=0, totalM=0, scale=0;
    char comma;
    if (!(ss >> target)) { ifs.close(); return -1.0; }
    ss >> comma;
    if (!(ss >> totalM)) { ifs.close(); return -1.0; }
    ss >> comma;
    if (!(ss >> scale)) { ifs.close(); return -1.0; }
    ifs.close();
    return scale;
}

void peakEnergy_converted_vs_beam(const char* muCSV = "/home/rudradeb/sim/output/muon_mip_calibration.csv",
                                  const char* scaleFile = "/home/rudradeb/sim/output/electron_scale.txt")
{
    gStyle->SetOptStat(0);

    // Edit these if your electron files/energies differ
    std::vector<int> energies = {10,20,50,100};
    std::map<int,std::string> fileOf = {
        {10,  "/home/rudradeb/sim/output/output_thread-2.root"},
        {20,  "/home/rudradeb/sim/output/output_thread-3.root"},
        {50,  "/home/rudradeb/sim/output/output_thread-4.root"},
        {100, "/home/rudradeb/sim/output/output_thread-5.root"}
    };

    // number of layers in your detector
    const int nLayers = 24;

    // read muon per-cell means
    double overallMean = 0.0;
    auto muMean = readMuonCSV(muCSV, overallMean);
    if (muMean.empty()) {
        printf("ERROR: muon calibration empty. Provide a valid muCSV (muon_mip_calibration.csv).\n");
        return;
    }
    printf("Read muon calibration: overallMean = %.6e (units of fEdep)\n", overallMean);

    // read global scale (GeV per MIP)
    double scale = readScale(scaleFile);
    if (scale <= 0) {
        printf("ERROR: cannot read global scale from %s. Provide electron_scale.txt with scale.\n", scaleFile);
        return;
    }
    printf("Read global scale: 1 MIP = %.12e GeV\n", scale);

    // output CSV
    std::string outCSV = "/home/rudradeb/sim/output/peak_energy_converted.csv";
    std::ofstream ofs(outCSV);
    if (!ofs.is_open()) { printf("ERROR: cannot write %s\n", outCSV.c_str()); return; }
    ofs << "EnergyGeV,nEvents,peakLayer,peakMeanEnergy_GeV\n";

    // containers for plotting
    std::vector<double> xE;         // beam energies
    std::vector<double> yPeakEgev;  // peak mean energy (GeV)
    std::vector<double> yPeakLayer; // peak layer index

    // loop energies
    for (int E : energies) {
        auto it = fileOf.find(E);
        if (it == fileOf.end()) { printf("No file specified for %d GeV — skipping\n", E); continue; }
        std::string fname = it->second;
        TFile f(fname.c_str());
        if (f.IsZombie()) { Warning("peakEnergy_converted_vs_beam","Cannot open %s — skipping", fname.c_str()); continue; }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) { Warning("peakEnergy_converted_vs_beam","No Hits tree in %s — skipping", fname.c_str()); f.Close(); continue; }

        // branches robustly
        Int_t detID = 0;
        Double_t edep = 0;
        Long64_t evID = -1;
        bool hasEventBranch = false;
        if (t->GetListOfBranches()->FindObject("fDetectorID")) t->SetBranchAddress("fDetectorID",&detID);
        else if (t->GetListOfBranches()->FindObject("detID")) t->SetBranchAddress("detID",&detID);
        else { Warning("peakEnergy_converted_vs_beam","No detector ID branch in %s", fname.c_str()); f.Close(); continue; }

        if (t->GetListOfBranches()->FindObject("fEdep")) t->SetBranchAddress("fEdep",&edep);
        else if (t->GetListOfBranches()->FindObject("Edep")) t->SetBranchAddress("Edep",&edep);
        else { Warning("peakEnergy_converted_vs_beam","No energy branch in %s", fname.c_str()); f.Close(); continue; }

        if (t->GetListOfBranches()->FindObject("fEvent")) { t->SetBranchAddress("fEvent",&evID); hasEventBranch=true; }
        else if (t->GetListOfBranches()->FindObject("event")) { t->SetBranchAddress("event",&evID); hasEventBranch=true; }
        else if (t->GetListOfBranches()->FindObject("evt")) { t->SetBranchAddress("evt",&evID); hasEventBranch=true; }

        // accumulate per-event per-layer energies (in GeV). Use map: eventID -> vector[layer]
        std::map<long long, std::vector<double>> evLayerGeV;
        Long64_t N = t->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            int L,cx,ry; decodeDet(detID, L, cx, ry);
            if (L<0 || L>=nLayers) continue;
            long long thisEv = hasEventBranch ? (long long)evID : (long long)i; // fallback: each entry unique
            if (!evLayerGeV.count(thisEv)) evLayerGeV[thisEv] = std::vector<double>(nLayers, 0.0);

            // convert edep -> number of MIPs using muon per-cell mean MIP:
            double meanMip = overallMean;
            auto itm = muMean.find({cx,ry});
            if (itm != muMean.end() && itm->second > 0) meanMip = itm->second;
            double nMips = (meanMip > 0) ? (edep / meanMip) : 0.0;
            double eg = nMips * scale; // energy in GeV for this hit
            evLayerGeV[thisEv][L] += eg;
        }
        f.Close();

        // compute mean energy per layer across events (GeV)
        int nEvents = evLayerGeV.size();
        if (nEvents == 0) { Warning("peakEnergy_converted_vs_beam","No events for %d GeV", E); continue; }
        std::vector<double> meanLayerGeV(nLayers, 0.0);
        for (auto &p : evLayerGeV) {
            const std::vector<double> &v = p.second;
            for (int L=0; L<nLayers; ++L) meanLayerGeV[L] += v[L];
        }
        for (int L=0; L<nLayers; ++L) meanLayerGeV[L] /= double(nEvents);

        // find peak layer & its mean energy (GeV)
        int peakL = std::distance(meanLayerGeV.begin(), std::max_element(meanLayerGeV.begin(), meanLayerGeV.end()));
        double peakMeanGeV = meanLayerGeV[peakL];

        printf("E=%d GeV: nEvents=%d peakLayer=%d  peakMeanEnergy = %.6e GeV\n", E, nEvents, peakL, peakMeanGeV);

        // store results
        xE.push_back(double(E));
        yPeakEgev.push_back(peakMeanGeV);
        yPeakLayer.push_back(double(peakL));

        ofs << E << "," << nEvents << "," << peakL << "," << std::setprecision(12) << peakMeanGeV << "\n";
    } // energies loop

    ofs.close();
    printf("Wrote CSV: %s\n", outCSV.c_str());

    if (xE.empty()) { printf("No points to plot.\n"); return; }

    // Draw plot: peak mean energy (GeV) vs beam energy (GeV), log-x recommended
    std::string outPNG = "/home/rudradeb/sim/output/peakEnergy_vs_beam.png";
    TCanvas *c = new TCanvas("c_peakE","Peak energy vs beam",900,600);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);

    // create frame with good ranges
    double xmin = *std::min_element(xE.begin(), xE.end());
    double xmax = *std::max_element(xE.begin(), xE.end());
    double ymin = *std::min_element(yPeakEgev.begin(), yPeakEgev.end());
    double ymax = *std::max_element(yPeakEgev.begin(), yPeakEgev.end());
    double xpad = 0.15*(xmax-xmin);
    double ypad = (ymax>ymin) ? 0.2*(ymax-ymin) : 0.1*ymax;
    TH1D *frame = new TH1D("frame","Peak mean energy in peak layer;Beam energy (GeV);Peak mean energy (GeV)", 10, xmin - xpad, xmax + xpad);
    frame->SetMinimum(std::max(0.0, ymin - ypad));
    frame->SetMaximum(ymax + ypad);
    frame->Draw();

    gPad->SetLogx(1); // beam energy on log x
    TGraph *g = new TGraph(xE.size());
    for (size_t i=0;i<xE.size();++i) g->SetPoint(i, xE[i], yPeakEgev[i]);
    g->SetMarkerStyle(21);
    g->SetMarkerSize(1.4);
    g->SetLineWidth(2);
    g->Draw("P SAME");

    // annotate points with layer index and value in GeV
    TLatex tx; tx.SetNDC(kFALSE); tx.SetTextSize(0.03);
    for (size_t i=0;i<xE.size();++i) {
        double px = xE[i];
        double py = yPeakEgev[i];
        // draw small label with layer and value
        TLatex t; t.SetTextSize(0.03);
        std::stringstream ss; ss << "L=" << int(yPeakLayer[i]) << "  " << std::fixed << std::setprecision(3) << yPeakEgev[i] << " GeV";
        // place text slightly above the point in data coordinates:
        t.DrawLatex(px*1.02, py*1.08, ss.str().c_str());
    }

    // title and save
    TLatex ttl; ttl.SetNDC(); ttl.SetTextSize(0.04);
    ttl.DrawLatex(0.30, 0.95, "Peak mean energy (GeV) vs beam energy");
    c->SaveAs(outPNG.c_str());
    printf("Saved plot: %s\n", outPNG.c_str());

    delete frame;
    delete g;
    delete c;
    printf("peakEnergy_converted_vs_beam: finished.\n");
}

