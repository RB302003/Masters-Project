// longitudinal_profile_muons.C
// Compute longitudinal (layer) profile from muon runs.
// Writes CSV and saves a PNG plot.
//
// Usage:
//   root -l longitudinal_profile_muons.C
//   longitudinal_profile_muons();                           // default args
//   longitudinal_profile_muons("/home/you/muon.root", 24);  // specify file & nLayers

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    // detID format assumed: layer*10000 + col*100 + row
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

void longitudinal_profile_muons(const char* muFile = "/home/rudradeb/sim/muonoutput/output_thread-5.root",
                                int nLayers = 24,
                                double timeMin = 0.0,
                                double timeMax = 1000.0,
                                const char* outCSV = "/home/rudradeb/sim/output/longitudinal_profile_muons.csv",
                                const char* outPNG = "/home/rudradeb/sim/output/longitudinal_profile_muons4.png")
{
    gStyle->SetOptStat(0);

    // Open file and tree
    TFile f(muFile);
    if (f.IsZombie()) {
        printf("ERROR: cannot open muon file: %s\n", muFile);
        return;
    }
    TTree *t = (TTree*)f.Get("Hits");
    if (!t) {
        printf("ERROR: Hits tree not found in %s\n", muFile);
        f.Close();
        return;
    }

    // Find branches robustly
    Int_t detID = 0;
    Double_t edep = 0;
    Int_t evID = -1;
    Double_t timeb = 0;
    bool hasEventBranch = false;
    bool hasTime = false;

    if (t->GetListOfBranches()->FindObject("fDetectorID")) t->SetBranchAddress("fDetectorID",&detID);
    else if (t->GetListOfBranches()->FindObject("detID")) t->SetBranchAddress("detID",&detID);
    else {
        printf("ERROR: no detector ID branch (fDetectorID/detID)\n"); f.Close(); return;
    }

    if (t->GetListOfBranches()->FindObject("fEdep")) t->SetBranchAddress("fEdep",&edep);
    else if (t->GetListOfBranches()->FindObject("Edep")) t->SetBranchAddress("Edep",&edep);
    else { printf("ERROR: no energy branch (fEdep/Edep)\n"); f.Close(); return; }

    if (t->GetListOfBranches()->FindObject("fEvent")) { t->SetBranchAddress("fEvent",&evID); hasEventBranch = true; }
    else if (t->GetListOfBranches()->FindObject("event")) { t->SetBranchAddress("event",&evID); hasEventBranch = true; }
    else if (t->GetListOfBranches()->FindObject("evt")) { t->SetBranchAddress("evt",&evID); hasEventBranch = true; }
    else {
        // Not fatal: we'll use one-entry-per-hit fallback if event id missing,
        // but better to have event-level accumulation; warn the user.
        printf("WARNING: no event ID branch found -> using entry index as event ID fallback.\n");
    }

    if (t->GetListOfBranches()->FindObject("time") ) { t->SetBranchAddress("time",&timeb); hasTime = true; }
    else if (t->GetListOfBranches()->FindObject("fTime") ) { t->SetBranchAddress("fTime",&timeb); hasTime = true; }

    // Map from eventID -> vector<double> layer sums (size nLayers)
    std::map<long long, std::vector<double>> evLayerSum;
    // Also keep track of events we encountered, to count events later
    std::map<long long, bool> evSeen;

    Long64_t N = t->GetEntries();
    for (Long64_t i=0; i<N; ++i) {
        t->GetEntry(i);

        if (edep <= 0) continue;
        if (hasTime && (timeb < timeMin || timeb > timeMax)) continue;

        long long thisEvent = 0;
        if (hasEventBranch) thisEvent = (long long)evID;
        else thisEvent = i; // fallback: treat each entry as its own event (not ideal but safe)

        int L, cx, ry;
        decodeDet(detID, L, cx, ry);
        if (L < 0 || L >= nLayers) continue;

        // ensure vector exists
        if (!evLayerSum.count(thisEvent)) evLayerSum[thisEvent] = std::vector<double>(nLayers, 0.0);

        evLayerSum[thisEvent][L] += edep;
        evSeen[thisEvent] = true;
    }

    // Now compute mean and RMS across events for each layer
    int nEvents = evLayerSum.size();
    if (nEvents == 0) { printf("No events collected after cuts.\n"); f.Close(); return; }
    printf("Events considered = %d   (time window = [%.3f, %.3f])\n", nEvents, timeMin, timeMax);

    std::vector<double> meanLayer(nLayers, 0.0);
    std::vector<double> rmsLayer(nLayers, 0.0);
    std::vector<int> nEventsWithSignal(nLayers, 0);

    // accumulate sums & sumsq per layer
    for (auto &p : evLayerSum) {
        const std::vector<double> &vec = p.second;
        for (int L=0; L<nLayers; ++L) {
            double v = vec[L];
            meanLayer[L] += v;
            rmsLayer[L] += v*v;
            if (v > 0) nEventsWithSignal[L] += 1;
        }
    }

    for (int L=0; L<nLayers; ++L) {
        meanLayer[L] /= double(nEvents);
        double meanSq = rmsLayer[L] / double(nEvents);
        double var = meanSq - meanLayer[L]*meanLayer[L];
        if (var < 0 && var > -1e-12) var = 0; // numerical safety
        rmsLayer[L] = (var > 0) ? sqrt(var) : 0.0;
    }

    // Write CSV
    {
        std::ofstream ofs(outCSV);
        ofs << "layer,meanE,rmsE,nEventsWithSignal\n";
        for (int L=0; L<nLayers; ++L) {
            ofs << L << "," << std::setprecision(12) << meanLayer[L] << "," << rmsLayer[L] << "," << nEventsWithSignal[L] << "\n";
        }
        ofs.close();
        printf("Wrote longitudinal profile CSV: %s\n", outCSV);
    }

    // Plot: TGraphErrors mean vs layer with RMS as error
    {
        TCanvas *c = new TCanvas("c_long", "Longitudinal profile (muons)", 900, 600);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);

        TGraphErrors *g = new TGraphErrors(nLayers);
        for (int L=0; L<nLayers; ++L) {
            double x = L + 1; // 1-based layer number on x-axis for nicer readability
            double y = meanLayer[L];
            double ey = rmsLayer[L] / sqrt( (nEvents>0) ? double(nEvents) : 1.0 ); // standard error of mean
            g->SetPoint(L, x, y);
            g->SetPointError(L, 0.0, ey);
        }

        // Autoscale Y limits with margin
        double ymax = 0;
        double ymin = 1e300;
        for (int L=0; L<nLayers; ++L) {
            if (meanLayer[L] > ymax) ymax = meanLayer[L];
            if (meanLayer[L] > 0 && meanLayer[L] < ymin) ymin = meanLayer[L];
        }
        if (ymin > ymax/2) ymin = 0; // safety
        double ypad = (ymax>0) ? 0.15*ymax : 1.0;

        // Create frame histogram for axes
        double xmin = 0.5, xmax = double(nLayers) + 0.5;
        TH1D *frame = new TH1D("frame_long", Form("Longitudinal profile (muons);Layer;Mean energy per event (fEdep units)"), nLayers, xmin, xmax);
        frame->SetMinimum(std::max(0.0, ymin - ypad));
        frame->SetMaximum(ymax + ypad);
        frame->Draw();

        // Draw graph
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetLineWidth(2);
        g->Draw("P same");

        // Add text summary
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.DrawLatex(0.65, 0.88, Form("Events = %d", nEvents));
        tx.DrawLatex(0.65, 0.82, Form("Time window = [%.3f, %.3f]", timeMin, timeMax));

        // Save PNG
        c->SaveAs(outPNG);
        printf("Saved longitudinal plot PNG: %s\n", outPNG);

        // clean
        delete g;
        delete frame;
        delete c;
    }

    f.Close();
    printf("longitudinal_profile_muons: finished.\n");
}

