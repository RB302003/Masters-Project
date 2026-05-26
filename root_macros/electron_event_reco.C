// electron_event_reco_validated.C
// Event-level reconstruction for electron runs using muon calibration + global scale.
// Adds an automatic validation step that checks the scale by reconstructing
// the average energy of the 10 GeV run and warning if the scale is inconsistent.
//
// Usage:
//   root -l electron_event_reco_validated.C
//   electron_event_reco_validated(); // uses defaults
//   electron_event_reco_validated(muCSV, scaleFile); // custom paths

#include <map>
#include <vector>
#include <string>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

// Read muon CSV: col,row,meanE,calibFactor,countHits
static std::map<std::pair<int,int>, double> readMuonCSV(const char* muCSV, double &overallMean) {
    std::map<std::pair<int,int>, double> muMean;
    overallMean = 0.0;
    std::ifstream ifs(muCSV);
    if (!ifs.is_open()) {
        printf("WARNING: cannot open muon CSV '%s' — muMean will be empty.\n", muCSV);
        return muMean;
    }
    std::string line;
    if (!std::getline(ifs, line)) { ifs.close(); return muMean; } // header
    double sum = 0.0; int count = 0;
    while (std::getline(ifs, line)) {
        if (line.size() == 0) continue;
        std::stringstream ss(line);
        int col, row; double meanE = 0.0; double calib = 0.0; long long cnt = 0;
        char comma;
        if (!(ss >> col)) continue;
        ss >> comma;
        if (!(ss >> row)) continue;
        ss >> comma;
        if (!(ss >> meanE)) meanE = 0.0;
        muMean[{col,row}] = meanE;
        if (meanE > 0) { sum += meanE; ++count; }
    }
    ifs.close();
    overallMean = (count>0) ? (sum/double(count)) : 0.0;
    return muMean;
}

// Read scale file format:
// targetEnergyGeV, N_events, total_MIPs_sum, scale_GeV_per_MIP
static double readScale(const char* scaleFile, double &scale_out) {
    scale_out = -1.0;
    std::ifstream ifs(scaleFile);
    if (!ifs.is_open()) {
        printf("WARNING: cannot open scale file '%s'\n", scaleFile);
        return -1.0;
    }
    std::string line;
    // skip header if present
    if (!std::getline(ifs, line)) { ifs.close(); return -1.0; }
    // try second line that contains values
    if (!std::getline(ifs, line)) { ifs.close(); return -1.0; }
    std::stringstream ss(line);
    double target=0.0; double nEvents=0.0; double totalM=0.0; double scale=0.0;
    char comma;
    if (!(ss >> target)) { ifs.close(); return -1.0; }
    ss >> comma;
    if (!(ss >> nEvents)) { ifs.close(); return -1.0; }
    ss >> comma;
    if (!(ss >> totalM)) { ifs.close(); return -1.0; }
    ss >> comma;
    if (!(ss >> scale)) { ifs.close(); return -1.0; }
    ifs.close();
    scale_out = scale;
    return scale;
}

// Main macro
void electron_event_reco(const char* muCSV = "/tmp/muon_mip_calibration.csv",
                                   const char* scaleFile = "/home/rudradeb/sim/output/electron_scale.txt")
{
    // Map of electron files (edit paths if needed)
    std::map<int,std::string> fileOf = {
        {1,  "/home/rudradeb/sim/finaloutput/output_thread-1.root"},
        {2,  "/home/rudradeb/sim/finaloutput/output_thread-2.root"},
        {5,  "/home/rudradeb/sim/finaloutput/output_thread-3.root"},
        {10, "/home/rudradeb/sim/finaloutput/output_thread-4.root"}
    };
    std::vector<int> energies = {1,2,5,10};
    const int nSide = 7;

    // 1) Read muon calibration
    double overallMean = 0.0;
    auto muMean = readMuonCSV(muCSV, overallMean);
    if (muMean.empty()) {
        printf("ERROR: muon calibration map empty. Run muon_calibrate and provide correct CSV.\n");
        return;
    }
    printf("Muon calibration: overall mean MIP = %.6e (fallback if a cell missing)\n", overallMean);

    // 2) Read scale and validate existence
    double scale = -1.0;
    if (readScale(scaleFile, scale) <= 0.0) {
        printf("ERROR: failed to read valid scale from '%s'. Aborting.\n", scaleFile);
        return;
    }
    printf("Global scale read: 1 MIP = %.12e GeV\n", scale);

    // ===============================
    // SCALE VALIDATION STEP (uses 10 GeV run)
    // ===============================
    {
        auto it10 = fileOf.find(10);
        if (it10 != fileOf.end()) {
            std::string valFile = it10->second;
            TFile vf(valFile.c_str());
            if (!vf.IsZombie()) {
                TTree *vt = (TTree*)vf.Get("Hits");
                if (vt) {
                    Int_t detID_val = 0, evID_val = 0;
                    Double_t edep_val = 0;
                    bool hasDet = false, hasEdep = false, hasEv = false;
                    if (vt->GetListOfBranches()->FindObject("fDetectorID")) { vt->SetBranchAddress("fDetectorID",&detID_val); hasDet=true; }
                    else if (vt->GetListOfBranches()->FindObject("detID")) { vt->SetBranchAddress("detID",&detID_val); hasDet=true; }
                    if (vt->GetListOfBranches()->FindObject("fEdep")) { vt->SetBranchAddress("fEdep",&edep_val); hasEdep=true; }
                    else if (vt->GetListOfBranches()->FindObject("Edep")) { vt->SetBranchAddress("Edep",&edep_val); hasEdep=true; }
                    if (vt->GetListOfBranches()->FindObject("fEvent")) { vt->SetBranchAddress("fEvent",&evID_val); hasEv=true; }
                    else if (vt->GetListOfBranches()->FindObject("event")) { vt->SetBranchAddress("event",&evID_val); hasEv=true; }
                    else if (vt->GetListOfBranches()->FindObject("evt")) { vt->SetBranchAddress("evt",&evID_val); hasEv=true; }

                    if (!hasDet || !hasEdep) {
                        printf("Validation: missing fDetectorID or fEdep branch in %s — skipping validation.\n", valFile.c_str());
                    } else {
                        std::map<int,double> mipsPerEvent;
                        Long64_t Nval = vt->GetEntries();
                        for (Long64_t i=0;i<Nval;++i) {
                            vt->GetEntry(i);
                            if (edep_val <= 0) continue;
                            int L,cx,ry; decodeDet(detID_val, L, cx, ry);
                            if (cx<0 || cx>=nSide || ry<0 || ry>=nSide) continue;
                            double meanMIP = overallMean;
                            auto it = muMean.find({cx,ry});
                            if (it != muMean.end() && it->second > 0) meanMIP = it->second;
                            double nM = (meanMIP>0) ? (edep_val / meanMIP) : 0.0;
                            int ev = hasEv ? evID_val : int(i); // fallback: unique per-entry id
                            mipsPerEvent[ev] += nM;
                        }
                        int Nevents = mipsPerEvent.size();
                        double sumM = 0.0;
                        for (auto &p : mipsPerEvent) sumM += p.second;
                        double avgMIP = (Nevents>0) ? (sumM / double(Nevents)) : 0.0;
                        double expectedE = avgMIP * scale;
                        printf("\n=== SCALE VALIDATION (10 GeV run) ===\n");
                        printf("Events found = %d\n", Nevents);
                        printf("Average MIPs/event (10 GeV run): %.6e\n", avgMIP);
                        printf("Reconstructed avg energy = avgMIP * scale = %.6f GeV\n", expectedE);
                        printf("Target beam energy = 10 GeV\n");
                        double ratio = (10.0>0) ? (expectedE / 10.0) : 0.0;
                        if (ratio < 0.8 || ratio > 1.2) {
                            printf("WARNING: Scale seems inconsistent! Reconstructed is %.3f × the beam energy.\n", ratio);
                        } else {
                            printf("Scale OK (within ±20%%).\n");
                        }
                        printf("=====================================\n\n");
                    }
                } else {
                    printf("Validation: Hits tree missing in %s — skipping validation.\n", valFile.c_str());
                }
            } else {
                printf("Validation: cannot open validation file %s — skipping validation.\n", valFile.c_str());
            }
            vf.Close();
        } else {
            printf("Validation: no 10 GeV filename configured — skipping validation.\n");
        }
    }
    // =============================== end validation

    // Prepare summary CSV
    const char *summaryCSV = "/home/rudradeb/sim/output/energy_summary.csv";
    FILE *sf = fopen(summaryCSV, "w");
    if (!sf) { printf("ERROR: cannot write summary CSV %s\n", summaryCSV); return; }
    fprintf(sf, "EnergyGeV, Nevents, MPV_fit_GeV, sigma_fit_GeV, mean_GeV, rms_GeV, resolution_fit, resolution_rms\n");

    // loop energies
    for (int E : energies) {
        auto itf = fileOf.find(E);
        if (itf == fileOf.end()) { printf("No file specified for %d GeV — skipping\n", E); continue; }
        std::string fname = itf->second;
        TFile f(fname.c_str());
        if (f.IsZombie()) { Warning("electron_event_reco_validated","Cannot open %s — skipping", fname.c_str()); continue; }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) { Warning("electron_event_reco_validated","No Hits tree in %s — skipping", fname.c_str()); f.Close(); continue; }

        // Branch discovery
        Int_t detID = 0;
        Double_t edep = 0;
        Int_t evID = 0;
        bool hasDet=false, hasEdep=false, hasEv=false;
        if (t->GetListOfBranches()->FindObject("fDetectorID")) { t->SetBranchAddress("fDetectorID",&detID); hasDet=true; }
        else if (t->GetListOfBranches()->FindObject("detID")) { t->SetBranchAddress("detID",&detID); hasDet=true; }
        if (t->GetListOfBranches()->FindObject("fEdep")) { t->SetBranchAddress("fEdep",&edep); hasEdep=true; }
        else if (t->GetListOfBranches()->FindObject("Edep")) { t->SetBranchAddress("Edep",&edep); hasEdep=true; }
        if (t->GetListOfBranches()->FindObject("fEvent")) { t->SetBranchAddress("fEvent",&evID); hasEv=true; }
        else if (t->GetListOfBranches()->FindObject("event")) { t->SetBranchAddress("event",&evID); hasEv=true; }
        else if (t->GetListOfBranches()->FindObject("evt")) { t->SetBranchAddress("evt",&evID); hasEv=true; }

        if (!hasDet || !hasEdep) { Warning("electron_event_reco_validated","Missing branches in %s — skipping", fname.c_str()); f.Close(); continue; }
        if (!hasEv) {
            printf("WARNING: no event ID branch found in %s — treating each hit as unique event (not ideal).\n", fname.c_str());
        }

        // accumulate per-event totals using a map eventID -> pair(totalE, totalMips)
        std::map<int, std::pair<double,double>> evMap; // evID -> (totalE_units, totalMips)
        Long64_t N = t->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            int L, cx, ry; decodeDet(detID, L, cx, ry);
            if (cx<0 || cx>=nSide || ry<0 || ry>=nSide) continue;
            double meanMIP = overallMean;
            auto it = muMean.find({cx,ry});
            if (it != muMean.end() && it->second > 0) meanMIP = it->second;
            double nMips_hit = (meanMIP > 0) ? (edep / meanMIP) : 0.0;
            int thisEvent = hasEv ? evID : int(i); // fallback: unique id per entry
            evMap[thisEvent].first += edep;
            evMap[thisEvent].second += nMips_hit;
        }
        f.Close();

        // collect per-event vectors
        std::vector<double> vecEgev;
        std::vector<double> vecMips;
        vecEgev.reserve(evMap.size());
        vecMips.reserve(evMap.size());
        // CSV for per-event
        char outEvents[256];
        snprintf(outEvents, sizeof(outEvents), "/home/rudradeb/sim/output/events_E%dGeV.csv", E);
        FILE *of = fopen(outEvents, "w");
        if (!of) { printf("ERROR: cannot write %s\n", outEvents); continue; }
        fprintf(of, "eventID,totalE_units,totalMIPs,energy_GeV\n");

        for (auto &p : evMap) {
            int id = p.first;
            double totE = p.second.first;
            double totM = p.second.second;
            double egev = totM * scale;
            vecEgev.push_back(egev);
            vecMips.push_back(totM);
            fprintf(of, "%d,%.12e,%.12e,%.12e\n", id, totE, totM, egev);
        }
        fclose(of);
        printf("Wrote per-event CSV: %s  (Nevents=%zu)\n", outEvents, vecEgev.size());

        if (vecEgev.empty()) { printf("No events found for %d GeV — skipping histogram/fit.\n", E); continue; }

        // Build histogram: choose binning around beam energy
        double histMin = 0.0;
        double histMax = 1.6 * double(E);
        int nBins = 200;
        TH1D *h = new TH1D(Form("h_E_%d",E), Form("Reconstructed energy - %d GeV;E (GeV);Entries", E), nBins, histMin, histMax);
        for (double val : vecEgev) h->Fill(val);
        double mean = h->GetMean();
        double rms = h->GetRMS();

        // Fit Landau around peak
        int maxBin = h->GetMaximumBin();
        double peakCenter = h->GetBinCenter(maxBin);
        double amp = h->GetBinContent(maxBin);
        double initMPV = peakCenter;
        double initSigma = rms > 0 ? rms/2.0 : 0.1*E;
        double fitMin = std::max(histMin, initMPV - 0.8*E);
        double fitMax = std::min(histMax, initMPV + 0.8*E);
        TF1 *fland = new TF1(Form("fland_%d",E), "landau", fitMin, fitMax);
        fland->SetParameters(amp, initMPV, initSigma);
        int fitStatus = h->Fit(fland, "RQ");
        double mpv = fland->GetParameter(1);
        double sigma_fit = fland->GetParameter(2);
        double res_fit = (mpv != 0) ? (sigma_fit / mpv) : 0.0;
        double res_rms = (mean != 0) ? (rms / mean) : 0.0;

        // Save histogram PNG
        TCanvas *c = new TCanvas(Form("c_E_%d",E), Form("Energy %d GeV",E), 900, 700);
        h->SetStats(0);
        h->Draw();
        fland->SetLineColor(kRed);
        fland->SetLineWidth(2);
        fland->Draw("same");
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
        tx.DrawLatex(0.62, 0.80, Form("MPV_{fit}=%.4f GeV", mpv));
        tx.DrawLatex(0.62, 0.74, Form("#sigma_{fit}=%.4f GeV", sigma_fit));
        tx.DrawLatex(0.62, 0.68, Form("mean=%.4f GeV  rms=%.4f GeV", mean, rms));
        tx.DrawLatex(0.62, 0.62, Form("res_{fit}=%.3f  res_{rms}=%.3f", res_fit, res_rms));
        char outpng[256];
        snprintf(outpng, sizeof(outpng), "/home/rudradeb/sim/output/energy_hist_E%dGeV.png", E);
        c->SaveAs(outpng);
        printf("Saved histogram+fit: %s\n", outpng);

        // Append to summary CSV
        fprintf(sf, "%d, %zu, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e\n",
                E, vecEgev.size(), mpv, sigma_fit, mean, rms, res_fit, res_rms);

        // cleanup
        delete h;
        delete fland;
        delete c;
    } // energies loop

    fclose(sf);
    printf("Summary CSV written to /home/rudradeb/sim/output/energy_summary.csv\n");
    printf("electron_event_reco_validated finished.\n");
}

