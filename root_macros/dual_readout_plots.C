/*
 * dual_readout_plots.C
 *
 * Generates the 4 dual-readout correlation plots from the PbWO4 homogeneous
 * calorimeter simulation output. Matches the style of the reference image:
 *   - S vs E_dep    (top-left)
 *   - C vs E_dep    (top-right)
 *   - C vs S        (bottom-left)
 *   - S+C vs E_dep  (bottom-right)
 *
 * Usage:
 *   root -l -b -q 'dual_readout_plots.C("output/output_thread0.root")'
 *
 * Output:
 *   dual_readout_plots.jpg  — 4-panel plot matching the reference image
 */

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TROOT.h"
#include <iostream>

void dual_readout_plots(const char* filename = "output/output_thread-2.root")
{
    /* ------------------------------------------------------------------ */
    /* 1. Style — match the clean ROOT default look from the reference     */
    /* ------------------------------------------------------------------ */
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1111);           // show entries, mean x/y, std dev x/y
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetPalette(kBird);          // default ROOT blue-green-yellow palette
    gStyle->SetNumberContours(64);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.16);   // room for COLZ colour bar
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.09);

    /* ------------------------------------------------------------------ */
    /* 2. Open ROOT file and load the Hits TTree                           */
    /* ------------------------------------------------------------------ */
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
        return;
    }

    TTree* hits = nullptr;
    f->GetObject("Hits", hits);
    if (!hits) {
        std::cerr << "ERROR: TTree 'Hits' not found in " << filename << std::endl;
        f->Close();
        return;
    }

    Long64_t nEntries = hits->GetEntries();
    std::cout << "Loaded TTree 'Hits' with " << nEntries << " entries." << std::endl;

    /* ------------------------------------------------------------------ */
    /* 3. Determine axis ranges from the data                              */
    /* ------------------------------------------------------------------ */
    Double_t eMin  = hits->GetMinimum("fEdep");
    Double_t eMax  = hits->GetMaximum("fEdep");
    Double_t sMin  = hits->GetMinimum("fScint");
    Double_t sMax  = hits->GetMaximum("fScint");
    Double_t cMin  = hits->GetMinimum("fCher");
    Double_t cMax  = hits->GetMaximum("fCher");
    Double_t scMax = sMax + cMax;

    /* ------------------------------------------------------------------ */
    /* 4. Book 2D histograms (100x100 bins, COLZ style)                   */
    /* ------------------------------------------------------------------ */
    // S vs E_dep
    TH2D* hE_S = new TH2D("hE_S",
        "S vs E_{dep};"
        "E_{dep} [MeV];"
        "S [a.u.]",
        100, eMin, eMax * 1.05,
        100, sMin, sMax * 1.05);

    // C vs E_dep
    TH2D* hE_C = new TH2D("hE_C",
        "C vs E_{dep};"
        "E_{dep} [MeV];"
        "C [a.u.]",
        100, eMin, eMax * 1.05,
        100, cMin, cMax * 1.05);

    // C vs S
    TH2D* hSC = new TH2D("hSC",
        "C vs S;"
        "S [a.u.];"
        "C [a.u.]",
        100, sMin, sMax * 1.05,
        100, cMin, cMax * 1.05);

    // S+C vs E_dep
    TH2D* hE_sum = new TH2D("hE_sum",
        "S+C vs E_{dep};"
        "E_{dep} [MeV];"
        "S+C [a.u.]",
        100, eMin, eMax * 1.05,
        100, 0,   scMax * 1.05);

    /* ------------------------------------------------------------------ */
    /* 5. Fill histograms by looping over the TTree                        */
    /* ------------------------------------------------------------------ */
    Double_t edep  = 0;
    Double_t scint = 0;
    Double_t cher  = 0;

    hits->SetBranchAddress("fEdep",  &edep);
    hits->SetBranchAddress("fScint", &scint);
    hits->SetBranchAddress("fCher",  &cher);

    for (Long64_t i = 0; i < nEntries; ++i) {
        hits->GetEntry(i);
        if (edep <= 0) continue;
        hE_S  ->Fill(edep, scint);
        hE_C  ->Fill(edep, cher);
        hSC   ->Fill(scint, cher);
        hE_sum->Fill(edep, scint + cher);
    }

    /* ------------------------------------------------------------------ */
    /* 6. Draw on a 2x2 canvas and save as JPG                            */
    /* ------------------------------------------------------------------ */
    TCanvas* c = new TCanvas("dual_readout", "Dual-Readout Correlation Plots",
                              1400, 1200);
    c->Divide(2, 2, 0.001, 0.001);

    // Top-left: S vs E_dep
    c->cd(1);
    gPad->SetLogz(0);
    hE_S->SetContour(64);
    hE_S->GetXaxis()->SetTitleOffset(1.2);
    hE_S->GetYaxis()->SetTitleOffset(1.4);
    hE_S->Draw("COLZ");

    // Top-right: C vs E_dep
    c->cd(2);
    gPad->SetLogz(0);
    hE_C->SetContour(64);
    hE_C->GetXaxis()->SetTitleOffset(1.2);
    hE_C->GetYaxis()->SetTitleOffset(1.4);
    hE_C->Draw("COLZ");

    // Bottom-left: C vs S
    c->cd(3);
    gPad->SetLogz(0);
    hSC->SetContour(64);
    hSC->GetXaxis()->SetTitleOffset(1.2);
    hSC->GetYaxis()->SetTitleOffset(1.4);
    hSC->Draw("COLZ");

    // Bottom-right: S+C vs E_dep
    c->cd(4);
    gPad->SetLogz(0);
    hE_sum->SetContour(64);
    hE_sum->GetXaxis()->SetTitleOffset(1.2);
    hE_sum->GetYaxis()->SetTitleOffset(1.4);
    hE_sum->Draw("COLZ");

    /* ------------------------------------------------------------------ */
    /* 7. Save outputs                                                     */
    /* ------------------------------------------------------------------ */
    c->Update();
    c->SaveAs("dual_readout_plots.jpg");
    std::cout << "\nSaved: dual_readout_plots.jpg" << std::endl;

    f->Close();
    delete f;
}
