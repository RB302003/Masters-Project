// make_muon_cal_plots_combined_withMeanEvsCalib.C
// Produce combined PNGs that overlay meanE per cell (log-y) and calib per cell,
// plus a MeanE vs Calib scatter (log-x) for all 4 muon runs.
// Writes per-run CSVs with mean and calib factor.
//
// Usage:
//   root -l make_muon_cal_plots_combined_withMeanEvsCalib.C
//   make_muon_cal_plots_combined_withMeanEvsCalib();

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

void make_muon_cal_plots(const char* baseDir = "/home/rudradeb/sim/muonoutput/",
                                  const char* outDir  = "/home/rudradeb/sim/output/",
                                  const char* filePattern = "output_thread-%d.root")
{
    gStyle->SetOptStat(0);
    const int nSide = 7;
    const int nCells = nSide * nSide;

    // list of 4 input files (adjust if you have different numbering)
    std::vector<std::string> files;
    for (int i = 1; i <= 4; ++i) {
        char buf[512];
        snprintf(buf, sizeof(buf), "%s", baseDir);
        std::string path = std::string(buf) + Form(filePattern, i);
        files.push_back(path);
    }

    // storage for per-file meanE and calib arrays
    std::vector<std::vector<double>> all_meanE;
    std::vector<std::vector<double>> all_calib;
    std::vector<std::string> labels;

    // color/marker palette
    int colors[4] = {kBlack, kRed+1, kBlue+1, kGreen+3};
    int markers[4] = {20, 21, 22, 23};

    // loop files, compute meanE per cell and calib (overallMean / cellMean)
    for (size_t fi=0; fi<files.size(); ++fi) {
        const std::string &muFile = files[fi];
        TFile f(muFile.c_str());
        if (f.IsZombie()) {
            printf("WARNING: cannot open muon file: %s  (skipping)\n", muFile.c_str());
            // push zeros so indices align
            all_meanE.emplace_back(nCells, 0.0);
            all_calib.emplace_back(nCells, 0.0);
            labels.push_back(Form("file %zu (missing)", fi+1));
            continue;
        }
        TTree *t = (TTree*)f.Get("Hits");
        if (!t) {
            printf("WARNING: Hits tree not found in %s (skipping)\n", muFile.c_str());
            all_meanE.emplace_back(nCells, 0.0);
            all_calib.emplace_back(nCells, 0.0);
            labels.push_back(Form("file %zu (no tree)", fi+1));
            f.Close();
            continue;
        }

        // branch bindings (robust)
        Int_t detID = 0;
        Double_t edep = 0;
        Double_t timeb = 0;
        bool hasTime = false;

        if (t->GetListOfBranches()->FindObject("fDetectorID")) t->SetBranchAddress("fDetectorID",&detID);
        else if (t->GetListOfBranches()->FindObject("detID")) t->SetBranchAddress("detID",&detID);
        else { printf("ERROR: no detector ID branch in %s\n", muFile.c_str()); f.Close(); return; }

        if (t->GetListOfBranches()->FindObject("fEdep")) t->SetBranchAddress("fEdep",&edep);
        else if (t->GetListOfBranches()->FindObject("Edep")) t->SetBranchAddress("Edep",&edep);
        else { printf("ERROR: no energy branch in %s\n", muFile.c_str()); f.Close(); return; }

        if (t->GetListOfBranches()->FindObject("time")) { t->SetBranchAddress("time",&timeb); hasTime=true; }
        else if (t->GetListOfBranches()->FindObject("fTime")) { t->SetBranchAddress("fTime",&timeb); hasTime=true; }

        // accumulate sums and counts per (col,row)
        std::vector<std::vector<double>> sumE(nSide, std::vector<double>(nSide,0.0));
        std::vector<std::vector<long long>> cnt(nSide, std::vector<long long>(nSide,0));

        Long64_t N = t->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            if (edep <= 0) continue;
            if (hasTime && (timeb < 0 || timeb > 1000)) continue;
            int L,cx,ry; decodeDet(detID, L, cx, ry);
            if (cx<0 || cx>=nSide || ry<0 || ry>=nSide) continue;
            sumE[cx][ry] += edep;
            cnt[cx][ry]  += 1;
        }

        // compute mean per cell and calib factor
        std::vector<double> meanE; meanE.reserve(nCells);
        std::vector<double> calib; calib.reserve(nCells);
        double sumMeans = 0.0; int nvalid = 0;
        for (int r=0;r<nSide;++r) {
            for (int c=0;c<nSide;++c) {
                double m = (cnt[c][r] > 0) ? (sumE[c][r] / double(cnt[c][r])) : 0.0;
                meanE.push_back(m);
                if (m > 0) { sumMeans += m; ++nvalid; }
            }
        }
        double overallMean = (nvalid>0) ? (sumMeans / double(nvalid)) : 0.0;
        for (size_t i=0;i<meanE.size();++i) calib.push_back(meanE[i] > 0 ? overallMean / meanE[i] : 0.0);

        // write CSV for this file
        std::string outCSV = std::string(outDir) + Form("muon_calibration_run%zu.csv", fi+1);
        {
            std::ofstream ofs(outCSV);
            ofs<<"col,row,meanE,calibFactor,countHits\n";
            for (int r=0, idx=0; r<nSide; ++r) {
                for (int c=0; c<nSide; ++c, ++idx) {
                    long long hits = cnt[c][r];
                    ofs<<c<<","<<r<<","<<std::setprecision(12)<<meanE[idx]<<","<<calib[idx]<<","<<hits<<"\n";
                }
            }
            ofs.close();
            printf("Wrote CSV: %s\n", outCSV.c_str());
        }

        // store
        all_meanE.push_back(meanE);
        all_calib.push_back(calib);
        labels.push_back(Form("run %zu", fi+1));

        f.Close();
    } // end files loop

    // ---------- Combined plot: meanE per cell for all runs (log-y) ----------
    {
        std::string outpng = std::string(outDir) + "meanE_per_cell_combined.png";
        TCanvas *C = new TCanvas("C_comb","MeanE per cell (combined)",1400,700);
        TH1D *frame = new TH1D("frame","MeanE per cell (combined)", nCells, -0.5, nCells-0.5);

        // determine ymin/ymax from all runs positive values
        double minPos = 1e99, maxVal = 0;
        for (size_t fi=0; fi<all_meanE.size(); ++fi) {
            for (double v : all_meanE[fi]) if (v>0) { minPos = std::min(minPos, v); maxVal = std::max(maxVal, v); }
        }
        if (minPos==1e99) minPos = 1e-6;
        double ymin = minPos*0.5;
        double ymax = (maxVal>0 ? maxVal*1.6 : 1.0);
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle("Cell index (row-major)");
        frame->GetYaxis()->SetTitle("Mean energy deposition (fEdep units)");
        frame->Draw();
        gPad->SetLogy(1);

        // create graphs for each run
        std::vector<TGraph*> graphs;
        for (size_t fi=0; fi<all_meanE.size(); ++fi) {
            TGraph *gr = new TGraph();
            gr->SetMarkerStyle(markers[fi % 4]);
            gr->SetMarkerColor(colors[fi % 4]);
            gr->SetLineColor(colors[fi % 4]);
            gr->SetMarkerSize(1.1);
            for (int i=0;i<nCells;i++) {
                double y = all_meanE[fi].size()==nCells ? (all_meanE[fi][i] > 0 ? all_meanE[fi][i] : ymin*0.75) : ymin*0.75;
                gr->SetPoint(gr->GetN(), i, y);
            }
            gr->Draw("P SAME");
            graphs.push_back(gr);
        }

        // x-axis labels (one set)
        // use label like "(c,r)"
        for (int i=0, idx=0;i<nSide;++i) {
            for (int j=0;j<nSide;++j, ++idx) {
                frame->GetXaxis()->SetBinLabel(idx+1, Form("(%d,%d)", j, i));
            }
        }
        frame->GetXaxis()->LabelsOption("v");

        // legend
        TLegend *leg = new TLegend(0.75,0.65,0.94,0.90);
        leg->SetFillColor(0); leg->SetBorderSize(0);
        for (size_t fi=0; fi<labels.size(); ++fi) {
            if (fi < graphs.size()) leg->AddEntry(graphs[fi], labels[fi].c_str(), "P");
            else leg->AddEntry((TObject*)0, labels[fi].c_str(), "");
        }
        leg->Draw();

        TLatex t; t.SetNDC(); t.SetTextSize(0.035);
        t.DrawLatex(0.02,0.96,"Mean energy per cell for 4 muon runs (overlay) - log y");

        C->SaveAs(outpng.c_str());
        printf("Saved %s\n", outpng.c_str());

        // cleanup
        delete frame;
        for (auto g: graphs) delete g;
        delete leg;
        delete C;
    }

    // ---------- Combined plot: calibration factor per cell for all runs (linear) ----------
    {
        std::string outpng = std::string(outDir) + "calib_per_cell_combined.png";
        TCanvas *C2 = new TCanvas("C2_comb","Calibration per cell (combined)",1400,700);
        TH1D *frame2 = new TH1D("frame2","Calib per cell (combined)", nCells, -0.5, nCells-0.5);

        // determine min/max calib across all runs
        double minC=1e99, maxC=0;
        for (size_t fi=0; fi<all_calib.size(); ++fi)
            for (double v: all_calib[fi]) if (v>0) { minC = std::min(minC,v); maxC = std::max(maxC,v); }
        if (minC==1e99) minC = 0.0;
        frame2->SetMinimum(std::max(0.0, minC*0.9 - 0.05));
        frame2->SetMaximum(maxC*1.2 + 0.05);
        frame2->GetXaxis()->SetTitle("Cell index (row-major)");
        frame2->GetYaxis()->SetTitle("Calibration factor");
        frame2->Draw();

        std::vector<TGraph*> graphs2;
        for (size_t fi=0; fi<all_calib.size(); ++fi) {
            TGraph *gr = new TGraph();
            gr->SetMarkerStyle(markers[fi % 4]);
            gr->SetMarkerColor(colors[fi % 4]);
            gr->SetLineColor(colors[fi % 4]);
            gr->SetMarkerSize(1.1);
            for (int i=0;i<nCells;i++) {
                double y = all_calib[fi].size()==nCells ? all_calib[fi][i] : 0.0;
                gr->SetPoint(gr->GetN(), i, y);
            }
            gr->Draw("P SAME");
            graphs2.push_back(gr);
        }

        for (int i=0, idx=0;i<nSide;++i) {
            for (int j=0;j<nSide;++j, ++idx) frame2->GetXaxis()->SetBinLabel(idx+1, Form("(%d,%d)", j, i));
        }
        frame2->GetXaxis()->LabelsOption("v");

        TLegend *leg2 = new TLegend(0.75,0.65,0.94,0.90);
        leg2->SetFillColor(0); leg2->SetBorderSize(0);
        for (size_t fi=0; fi<labels.size(); ++fi) {
            if (fi < graphs2.size()) leg2->AddEntry(graphs2[fi], labels[fi].c_str(), "P");
            else leg2->AddEntry((TObject*)0, labels[fi].c_str(), "");
        }
        leg2->Draw();

        TLatex t2; t2.SetNDC(); t2.SetTextSize(0.035);
        t2.DrawLatex(0.02,0.96,"Calibration factor per cell for 4 muon runs (overlay)");

        C2->SaveAs(outpng.c_str());
        printf("Saved %s\n", outpng.c_str());

        delete frame2;
        for (auto g: graphs2) delete g;
        delete leg2;
        delete C2;
    }

    // ---------- New: Combined MeanE vs Calib scatter (log-x) for all runs ----------
    {
        std::string outpng = std::string(outDir) + "meanE_vs_calib_combined.png";
        TCanvas *C3 = new TCanvas("C3_comb","MeanE vs Calib (combined)",1000,800);
        gPad->SetLogx(1);

        // collect points from all runs
        std::vector<double> xm_all, ym_all;
        std::vector<int> col_runIdx; // track run index for coloring
        for (size_t fi=0; fi<all_meanE.size(); ++fi) {
            for (size_t i=0;i<all_meanE[fi].size(); ++i) {
                double me = all_meanE[fi][i];
                double cb = all_calib[fi][i];
                if (me>0 && cb>0) {
                    xm_all.push_back(me);
                    ym_all.push_back(cb);
                    col_runIdx.push_back((int)fi);
                }
            }
        }

        if (xm_all.empty()) {
            printf("No positive meanE/calib pairs found — skipping MeanE vs Calib plot\n");
        } else {
            // axis ranges
            double xmin = *std::min_element(xm_all.begin(), xm_all.end());
            double xmax = *std::max_element(xm_all.begin(), xm_all.end());
            double ymin = *std::min_element(ym_all.begin(), ym_all.end());
            double ymax = *std::max_element(ym_all.begin(), ym_all.end());
            xmin *= 0.8; xmax *= 1.2;
            double yrmin = std::max(0.0, ymin*0.9 - 0.01);
            double yrmax = ymax*1.2 + 0.01;

            TH1D *frame3 = new TH1D("frame3","MeanE vs Calib", 10, xmin, xmax);
            frame3->SetMinimum(yrmin);
            frame3->SetMaximum(yrmax);
            frame3->GetXaxis()->SetTitle("Mean energy deposition (fEdep units) [log scale]");
            frame3->GetYaxis()->SetTitle("Calibration factor");
            frame3->Draw();

            // jittered scatter, colored by run
            std::mt19937_64 rng(1234567);
            std::normal_distribution<double> dist(0.0, 0.02);

            std::vector<TGraph*> scatters;
            for (int fi=0; fi<(int)files.size(); ++fi) {
                TGraph *g = new TGraph();
                g->SetMarkerStyle(markers[fi % 4]);
                g->SetMarkerColor(colors[fi % 4]);
                g->SetMarkerSize(1.1);
                scatters.push_back(g);
            }
            for (size_t i=0;i<xm_all.size();++i) {
                double jitter = std::pow(10.0, dist(rng));
                double xj = xm_all[i] * jitter;
                int rid = col_runIdx[i];
                scatters[rid]->SetPoint(scatters[rid]->GetN(), xj, ym_all[i]);
            }
            for (auto g: scatters) { g->Draw("P SAME"); }

            // binned mean line (log-bins)
            int nbins = std::min(12, (int)std::max(3, (int)xm_all.size()/3));
            std::vector<double> bxc, byc;
            double logmin = log10(xmin), logmax = log10(xmax);
            for (int b=0;b<nbins;b++) {
                double l0 = logmin + (logmax-logmin)*b/nbins;
                double l1 = logmin + (logmax-logmin)*(b+1)/nbins;
                double left = pow(10.0,l0), right = pow(10.0,l1);
                double sum = 0; int ct = 0;
                for (size_t i=0;i<xm_all.size();++i) if (xm_all[i]>=left && xm_all[i]<right) { sum += ym_all[i]; ++ct; }
                if (ct>0) {
                    bxc.push_back(pow(10.0,(l0+l1)/2.0));
                    byc.push_back(sum/ct);
                }
            }
            if (!bxc.empty()) {
                TGraph *gbin = new TGraph(bxc.size());
                for (size_t i=0;i<bxc.size();++i) gbin->SetPoint(i, bxc[i], byc[i]);
                gbin->SetLineColor(kBlack);
                gbin->SetLineWidth(2);
                gbin->SetMarkerStyle(21);
                gbin->Draw("LP SAME");
                delete gbin;
            }

            // Pearson correlation (log(x) vs y) and display
            double mx=0,my=0; int n = xm_all.size();
            std::vector<double> lx(n);
            for (int i=0;i<n;++i) { lx[i] = log(xm_all[i]); mx += lx[i]; my += ym_all[i]; }
            mx /= n; my /= n;
            double num=0, denx=0, deny=0;
            for (int i=0;i<n;++i) { double dx = lx[i]-mx; double dy = ym_all[i]-my; num += dx*dy; denx += dx*dx; deny += dy*dy; }
            double r = (denx>0 && deny>0) ? num / sqrt(denx*deny) : 0.0;
            TLatex tt; tt.SetNDC(); tt.SetTextSize(0.04);
            tt.DrawLatex(0.60,0.88,Form("Pearson r (log-x vs y) = %.3f", r));

            // legend for runs
            TLegend *leg = new TLegend(0.15,0.70,0.35,0.90);
            leg->SetFillColor(0); leg->SetBorderSize(0);
            for (size_t fi=0; fi<labels.size(); ++fi) {
                leg->AddEntry(scatters[fi], labels[fi].c_str(), "P");
            }
            leg->Draw();

            TLatex t; t.SetNDC(); t.SetTextSize(0.035);
            t.DrawLatex(0.02,0.96,"MeanE vs Calibration Factor (overlay runs)");

            C3->SaveAs(outpng.c_str());
            printf("Saved %s\n", outpng.c_str());

            // cleanup
            for (auto g: scatters) delete g;
            delete frame3;
            delete leg;
        }
        delete C3;
    }

    printf("Done: created combined overlays and CSVs in %s\n", outDir);
}

