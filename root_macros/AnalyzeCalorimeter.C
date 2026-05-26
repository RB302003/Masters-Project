#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace std;

// small helper to decode detID: layer = first 2 digits (detID/10000), col = middle, row = last
static inline void decodeDetID(int detID, int &layer, int &col, int &row) {
    layer = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int local = detID % 10000;
    col = local / 100;
    row = local % 100;
}

// small moving-average smoother (in-place) with window radius r
void movingAverageSmooth(vector<double> &v, int r=1) {
    if (r <= 0) return;
    int n = v.size();
    vector<double> tmp(n,0);
    for (int i=0;i<n;++i) {
        int a = max(0, i-r);
        int b = min(n-1, i+r);
        double s=0; int cnt=0;
        for (int j=a;j<=b;++j) { s += v[j]; ++cnt; }
        tmp[i] = (cnt>0 ? s/cnt : 0);
    }
    v.swap(tmp);
}

void AnalyzeCalorimeter() {
    gStyle->SetOptStat(0);

    // --- File list: do NOT change filenames (user requested) ---
    vector<string> files = {
        "/home/rudradeb/sim/finaloutput/output_thread-1.root",
        "/home/rudradeb/sim/finaloutput/output_thread-2.root",
        "/home/rudradeb/sim/finaloutput/output_thread-3.root",
        "/home/rudradeb/sim/finaloutput/output_thread-4.root"
    };

    // number of layers and lateral geometry (7x7)
    const int nLayers = 24;
    const int nSide = 7;
    const int colMin = 0, colMax = 6, rowMin = 0, rowMax = 6;

    // loop files one-by-one and produce outputs per file
    for (size_t fi=0; fi<files.size(); ++fi) {
        const string fname = files[fi];
        TFile *f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            cout << "⚠️ Could not open file: " << fname << "  (skipping)\n";
            if (f) { f->Close(); delete f; }
            continue;
        }
        cout << "📂 Processing file: " << fname << "\n";

        TTree *t = (TTree*) f->Get("Hits");
        if (!t) {
            cout << "⚠️ No 'Hits' tree in " << fname << "  (skipping)\n";
            f->Close(); delete f;
            continue;
        }

        // Branch variables (match your tree)
        Int_t fEvent = -1;
        Int_t fDetectorID = -1;
        Double_t fEdep = 0.0;
        Int_t fTrackID = -1;

        t->SetBranchAddress("fEvent", &fEvent);
        t->SetBranchAddress("fDetectorID", &fDetectorID);
        t->SetBranchAddress("fEdep", &fEdep);
        // track ID is present in your tree
        t->SetBranchAddress("fTrackID", &fTrackID);

        // ---- Step A: collect per-event, per-layer energy sums ----
        // We'll build: map eventID -> array[24] of layer sums
        std::map<int, std::vector<double>> eventLayerSum;
        std::set<int> eventIDs; // to count total events
        Long64_t nentries = t->GetEntries();
        for (Long64_t ent=0; ent<nentries; ++ent) {
            t->GetEntry(ent);
            if (fEdep <= 0) continue;
            int layer,col,row;
            decodeDetID(fDetectorID, layer, col, row);
            if (layer < 0 || layer >= nLayers) continue;

            if (eventLayerSum.find(fEvent) == eventLayerSum.end()) {
                eventLayerSum[fEvent] = std::vector<double>(nLayers, 0.0);
            }
            eventLayerSum[fEvent][layer] += fEdep;
            eventIDs.insert(fEvent);
        }

        // Number of events recorded in the file (based on fEvent presence)
        const size_t nEvents = eventIDs.size();
        if (nEvents == 0) {
            cout << "⚠️ No events found in " << fname << " (no fEvent values). Skipping.\n";
            f->Close(); delete f;
            continue;
        }
        cout << "   -> found " << nEvents << " unique events (based on fEvent)\n";

        // ---- Step B: compute mean energy per layer (average over events) ----
        vector<double> meanEdepPerLayer(nLayers, 0.0);
        for (int L=0; L<nLayers; ++L) {
            double sumOverEvents = 0.0;
            for (auto &evkv : eventLayerSum) {
                sumOverEvents += evkv.second[L];
            }
            meanEdepPerLayer[L] = sumOverEvents / double(nEvents);
        }

        // Make a smoothed copy for plotting to improve "bell" appearance
        vector<double> smoothLayer = meanEdepPerLayer;
        movingAverageSmooth(smoothLayer, 1); // window radius 1 -> 3-point moving average
        // You can increase radius (e.g. 2) if needed, but don't over-smooth.

        // --- Plot longitudinal profile (mean per layer) ---
        TCanvas *cLong = new TCanvas(Form("cLong_%zu", fi), Form("Longitudinal_%zu", fi), 900,600);
        TH1D *hLong = new TH1D(Form("hLong_%zu",fi), "Longitudinal profile;Layer;Mean Edep per event (MeV)", nLayers, -0.5, nLayers-0.5);
        for (int L=0; L<nLayers; ++L) hLong->SetBinContent(L+1, meanEdepPerLayer[L]);
        hLong->SetLineWidth(2);
        hLong->SetMarkerStyle(20);
        hLong->Draw("HIST");
        // overlay smoothed curve
        TGraph *gSmooth = new TGraph(nLayers);
        for (int L=0; L<nLayers; ++L) gSmooth->SetPoint(L, L, smoothLayer[L]);
        gSmooth->SetLineWidth(2);
        gSmooth->SetLineStyle(2);
        gSmooth->Draw("L SAME");

        // small legend
        TLegend *leg = new TLegend(0.65,0.75,0.88,0.88);
        leg->AddEntry(hLong, "Mean per layer (raw)", "l");
        leg->AddEntry(gSmooth, "Smoothed (3-bin MA)", "l");
        leg->Draw();

        cLong->SaveAs(Form("longitudinal_profile_file%zu.png", fi));
        cout << "   -> saved longitudinal_profile_file" << fi << ".png\n";

        // ---- Step C: find peak layer (shower max) using the mean (or total) ----
        int peakLayer = 0;
        double maxVal = -1;
        for (int L=0; L<nLayers; ++L) {
            if (meanEdepPerLayer[L] > maxVal) { maxVal = meanEdepPerLayer[L]; peakLayer = L; }
        }
        cout << "   -> peak (shower-max) layer = " << peakLayer << "   (mean Edep = " << maxVal << ")\n";

        // ---- Step D: build 7x7 lateral heatmap at peakLayer (accumulate over all events/files) ----
        // Re-loop tree entries and sum edep per (col,row) only for peakLayer
        // We'll also keep per-track totals (trackID -> edep) for the peak layer
        vector<vector<double>> cellSum(nSide, vector<double>(nSide, 0.0));
        std::map<int,double> trackEdepMap; // trackID -> edep (for peakLayer)

        t->SetBranchAddress("fEvent", &fEvent);
        t->SetBranchAddress("fDetectorID", &fDetectorID);
        t->SetBranchAddress("fEdep", &fEdep);
        t->SetBranchAddress("fTrackID", &fTrackID);
        for (Long64_t ent=0; ent<nentries; ++ent) {
            t->GetEntry(ent);
            if (fEdep <= 0) continue;
            int L,col,row; decodeDetID(fDetectorID, L, col, row);
            if (L != peakLayer) continue;
            if (col<colMin||col>colMax||row<rowMin||row>rowMax) continue;
            cellSum[col][row] += fEdep;
            trackEdepMap[fTrackID] += fEdep;
        }

        // create TH2D heatmap (column vs row)
        TH2D *hHeat = new TH2D(Form("hHeat_%zu",fi),
                               Form("Lateral heatmap (layer %d) - file %zu;Column;Row", peakLayer, fi),
                               nSide, colMin-0.5, colMax+0.5,
                               nSide, rowMin-0.5, rowMax+0.5);
        double globalMax = 0.0;
        for (int c=0;c<nSide;++c) for (int r=0;r<nSide;++r) if (cellSum[c][r] > globalMax) globalMax = cellSum[c][r];
        for (int c=0;c<nSide;++c) {
            for (int r=0;r<nSide;++r) {
                double val = cellSum[c][r];
                hHeat->SetBinContent(hHeat->GetXaxis()->FindBin(c), hHeat->GetYaxis()->FindBin(r), val);
            }
        }

        TCanvas *cHeat = new TCanvas(Form("cHeat_%zu",fi), "Lateral heatmap", 700,600);
        hHeat->SetContour(50);
        hHeat->Draw("COLZ");
        // annotate global max
        TPaveText *pt = new TPaveText(0.15,0.88,0.6,0.98,"NDC");
        pt->AddText(Form("file: %s", fname.c_str()));
        pt->AddText(Form("peak layer = %d  (mean Edep = %.3g)", peakLayer, maxVal));
        pt->Draw();
        cHeat->SaveAs(Form("lateral_heatmap_file%zu.png", fi));
        cout << "   -> saved lateral_heatmap_file" << fi << ".png\n";

        // ---- Step E: write top contributing trackIDs for this peak layer ----
        // produce a sorted vector of track contributions
        vector<pair<int,double>> trackVec;
        for (auto &kv : trackEdepMap) trackVec.push_back(kv);
        sort(trackVec.begin(), trackVec.end(), [](const pair<int,double>& a, const pair<int,double>& b){
            return a.second > b.second;
        });
        // write top 10 to a text file
        ofstream ofs(Form("top_tracks_file%zu.txt", fi));
        ofs << "# Top track contributions (file: " << fname << ", peakLayer=" << peakLayer << ")\n";
        ofs << "# trackID\tEdep(MeV)\tfraction_of_layer_total\n";
        // compute layer total
        double layerTotal = 0.0; for (auto &kv : cellSum) for (auto v : kv) layerTotal += v;
        for (size_t k=0; k < min(size_t(10), trackVec.size()); ++k) {
            int tid = trackVec[k].first;
            double e = trackVec[k].second;
            double frac = (layerTotal > 0 ? (e / layerTotal) : 0.0);
            ofs << tid << "\t" << e << "\t" << frac << "\n";
        }
        ofs.close();
        cout << "   -> wrote top_tracks_file" << fi << ".txt (layer total = " << layerTotal << " MeV)\n";

        // cleanup and close file
        f->Close(); delete f;

        // keep canvas objects alive until user closes ROOT; if running batch, they are saved already
        // (we don't delete the histograms to allow user to inspect them interactively)
    } // end file loop

    cout << "✅ Done. Generated longitudinal and lateral plots and top-track summaries for each file.\n";
    cout << "Note: PDG / particle-type is NOT available in your Hits tree. To get particle species,\n"
         << "you must record PDG code at hit creation in the simulation and then this macro can\n"
         << "be extended to summarize by PDG (electron/gamma/muon/pion/etc.).\n";
}

