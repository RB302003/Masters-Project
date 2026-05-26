// analyze_optical.C
// ROOT macro for comprehensive optical property analysis of PbWO4 calorimeter
// Usage: root -l -b -q 'analyze_optical.C("output/output_thread-1.root")'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void analyze_optical(const char* filename = "output/output_thread-1.root") {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    TTree* tOpt = (TTree*)file->Get("Optical");
    TTree* tPh  = (TTree*)file->Get("OpticalPhotons");
    TTree* tHits = (TTree*)file->Get("Hits");

    if (!tOpt || !tPh) {
        std::cerr << "Error: Optical trees not found in file" << std::endl;
        return;
    }

    std::cout << "=== Optical Analysis ===" << std::endl;
    std::cout << "Optical tree entries: " << tOpt->GetEntries() << std::endl;
    std::cout << "OpticalPhotons entries: " << tPh->GetEntries() << std::endl;
    if (tHits) std::cout << "Hits tree entries: " << tHits->GetEntries() << std::endl;

    // Create directories for the JPG groups
    gSystem->Exec("mkdir -p plots_jpg/production");
    gSystem->Exec("mkdir -p plots_jpg/spectrum");
    gSystem->Exec("mkdir -p plots_jpg/detection");

    // Open multi-page PDF
    TCanvas* cPdf = new TCanvas("cPdf", "PDF Wrapper", 1200, 500);
    cPdf->Print("optical_results.pdf[");

    // ======================================================================
    // PLOT 1: Cherenkov vs Scintillation Photon Production (per event)
    // ======================================================================
    TCanvas* c1 = new TCanvas("c1", "Photon Production", 1200, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    TH1D* hCher = new TH1D("hCher", "", 100, 0, 0);
    TH1D* hScint = new TH1D("hScint", "", 100, 0, 0);
    tOpt->Draw("nCherenkov>>hCher", "", "goff");
    tOpt->Draw("nScintillation>>hScint", "", "goff");

    // Auto-range
    double maxCher = hCher->GetXaxis()->GetXmax();
    double maxScint = hScint->GetXaxis()->GetXmax();

    delete hCher; delete hScint;

    hCher = new TH1D("hCher", ";Number of Photons;Events", 50, 0, maxCher * 1.2);
    hScint = new TH1D("hScint", ";Number of Photons;Events", 50, 0, maxScint * 1.2);
    tOpt->Draw("nCherenkov>>hCher", "", "goff");
    tOpt->Draw("nScintillation>>hScint", "", "goff");

    hCher->SetLineColor(kBlue);
    hCher->SetFillColorAlpha(kBlue, 0.3);
    hScint->SetLineColor(kRed);
    hScint->SetFillColorAlpha(kRed, 0.3);

    double ymax = TMath::Max(hCher->GetMaximum(), hScint->GetMaximum()) * 1.3;
    hCher->SetMaximum(ymax);

    hCher->Draw("HIST");
    hScint->Draw("HIST SAME");

    TLegend* leg1 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg1->AddEntry(hCher, "Cherenkov", "f");
    leg1->AddEntry(hScint, "Scintillation", "f");
    leg1->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.92, "Optical Photon Production per Event");

    // Scatter plot: Cherenkov vs Scintillation
    c1->cd(2);
    tOpt->Draw("nScintillation:nCherenkov", "", "colz");

    latex.DrawLatex(0.15, 0.92, "Cherenkov vs Scintillation Correlation");

    c1->SaveAs("plots_jpg/production/optical_production.jpg");
    c1->Print("optical_results.pdf");

    // ======================================================================
    // PLOT 2: Detected Photon Wavelength Spectrum
    // ======================================================================
    TCanvas* c2 = new TCanvas("c2", "Wavelength Spectrum", 1200, 500);
    c2->Divide(2, 1);

    c2->cd(1);
    TH1D* hWlAll = new TH1D("hWlAll", ";Wavelength [nm];Detected Photons", 100, 200, 800);
    TH1D* hWlCher = new TH1D("hWlCher", ";Wavelength [nm];Detected Photons", 100, 200, 800);
    TH1D* hWlScint = new TH1D("hWlScint", ";Wavelength [nm];Detected Photons", 100, 200, 800);

    tPh->Draw("wavelength>>hWlAll", "", "goff");
    tPh->Draw("wavelength>>hWlCher", "creatorProcess==0", "goff");
    tPh->Draw("wavelength>>hWlScint", "creatorProcess==1", "goff");

    hWlAll->SetLineColor(kBlack);
    hWlAll->SetLineWidth(2);
    hWlCher->SetLineColor(kBlue);
    hWlCher->SetFillColorAlpha(kBlue, 0.2);
    hWlScint->SetLineColor(kRed);
    hWlScint->SetFillColorAlpha(kRed, 0.2);

    hWlAll->Draw("HIST");
    hWlCher->Draw("HIST SAME");
    hWlScint->Draw("HIST SAME");

    TLegend* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg2->AddEntry(hWlAll, "All Detected", "l");
    leg2->AddEntry(hWlCher, "Cherenkov", "f");
    leg2->AddEntry(hWlScint, "Scintillation", "f");
    leg2->Draw();

    latex.DrawLatex(0.15, 0.92, "Detected Photon Wavelength Spectrum");

    // Detection time distribution
    c2->cd(2);
    TH1D* hTime = new TH1D("hTime", ";Detection Time [ns];Detected Photons", 100, 0, 100);
    TH1D* hTimeCher = new TH1D("hTimeCher", ";Detection Time [ns];Detected Photons", 100, 0, 100);
    TH1D* hTimeScint = new TH1D("hTimeScint", ";Detection Time [ns];Detected Photons", 100, 0, 100);

    tPh->Draw("time>>hTime", "", "goff");
    tPh->Draw("time>>hTimeCher", "creatorProcess==0", "goff");
    tPh->Draw("time>>hTimeScint", "creatorProcess==1", "goff");

    hTime->SetLineColor(kBlack);
    hTime->SetLineWidth(2);
    hTimeCher->SetLineColor(kBlue);
    hTimeCher->SetFillColorAlpha(kBlue, 0.2);
    hTimeScint->SetLineColor(kRed);
    hTimeScint->SetFillColorAlpha(kRed, 0.2);

    hTime->Draw("HIST");
    hTimeCher->Draw("HIST SAME");
    hTimeScint->Draw("HIST SAME");

    TLegend* leg3 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg3->AddEntry(hTime, "All Detected", "l");
    leg3->AddEntry(hTimeCher, "Cherenkov", "f");
    leg3->AddEntry(hTimeScint, "Scintillation", "f");
    leg3->Draw();

    latex.DrawLatex(0.15, 0.92, "Detection Time Distribution");

    c2->SaveAs("plots_jpg/spectrum/optical_spectrum.jpg");
    c2->Print("optical_results.pdf");

    // ======================================================================
    // PLOT 3: Detection Efficiency & Spatial Distribution
    // ======================================================================
    TCanvas* c3 = new TCanvas("c3", "Detection", 1200, 500);
    c3->Divide(2, 1);

    // Detection efficiency per event
    c3->cd(1);
    Int_t nEv, nCh, nSc, nDet, nDCh, nDSc;
    tOpt->SetBranchAddress("fEvent", &nEv);
    tOpt->SetBranchAddress("nCherenkov", &nCh);
    tOpt->SetBranchAddress("nScintillation", &nSc);
    tOpt->SetBranchAddress("nDetected", &nDet);
    tOpt->SetBranchAddress("nDetCherenkov", &nDCh);
    tOpt->SetBranchAddress("nDetScintillation", &nDSc);

    TH1D* hEff = new TH1D("hEff", ";Detection Efficiency;Events", 50, 0, 0.5);
    TH1D* hEffCher = new TH1D("hEffCher", ";Detection Efficiency;Events", 50, 0, 0.5);
    TH1D* hEffScint = new TH1D("hEffScint", ";Detection Efficiency;Events", 50, 0, 0.5);

    for (Long64_t i = 0; i < tOpt->GetEntries(); i++) {
        tOpt->GetEntry(i);
        G4int nTotal = nCh + nSc;
        if (nTotal > 0) {
            hEff->Fill((double)nDet / nTotal);
        }
        if (nCh > 0) {
            hEffCher->Fill((double)nDCh / nCh);
        }
        if (nSc > 0) {
            hEffScint->Fill((double)nDSc / nSc);
        }
    }

    hEff->SetLineColor(kBlack);
    hEff->SetLineWidth(2);
    hEffCher->SetLineColor(kBlue);
    hEffCher->SetFillColorAlpha(kBlue, 0.2);
    hEffScint->SetLineColor(kRed);
    hEffScint->SetFillColorAlpha(kRed, 0.2);

    double effMax = TMath::Max(hEff->GetMaximum(),
                     TMath::Max(hEffCher->GetMaximum(), hEffScint->GetMaximum())) * 1.3;
    hEff->SetMaximum(effMax);
    hEff->Draw("HIST");
    hEffCher->Draw("HIST SAME");
    hEffScint->Draw("HIST SAME");

    TLegend* leg4 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg4->AddEntry(hEff, "Overall", "l");
    leg4->AddEntry(hEffCher, "Cherenkov", "f");
    leg4->AddEntry(hEffScint, "Scintillation", "f");
    leg4->Draw();

    latex.DrawLatex(0.15, 0.92, "Photon Detection Efficiency");

    // Spatial distribution of detected photons (rear face X-Y)
    c3->cd(2);
    TH2D* hXY = new TH2D("hXY", ";X [mm];Y [mm]", 50, -80, 80, 50, -80, 80);
    tPh->Draw("posY:posX>>hXY", "", "colz goff");
    hXY->Draw("colz");

    latex.DrawLatex(0.15, 0.92, "Detected Photon Position (Rear Face)");

    c3->SaveAs("plots_jpg/detection/optical_detection.jpg");
    c3->Print("optical_results.pdf");

    // ======================================================================
    // PLOT 4: Summary Statistics
    // ======================================================================
    TCanvas* c4 = new TCanvas("c4", "Summary", 800, 600);

    // Print summary statistics
    std::cout << "\n=== OPTICAL PROPERTIES SUMMARY ===" << std::endl;

    double meanCher = 0, meanScint = 0, meanDet = 0;
    double meanEffTotal = 0;
    int nEvents = tOpt->GetEntries();

    for (Long64_t i = 0; i < nEvents; i++) {
        tOpt->GetEntry(i);
        meanCher += nCh;
        meanScint += nSc;
        meanDet += nDet;
        int total = nCh + nSc;
        if (total > 0) meanEffTotal += (double)nDet / total;
    }
    if (nEvents > 0) {
        meanCher /= nEvents;
        meanScint /= nEvents;
        meanDet /= nEvents;
        meanEffTotal /= nEvents;
    }

    std::cout << "Total events:              " << nEvents << std::endl;
    std::cout << "Mean Cherenkov/event:       " << meanCher << std::endl;
    std::cout << "Mean Scintillation/event:   " << meanScint << std::endl;
    std::cout << "Mean Detected/event:        " << meanDet << std::endl;
    std::cout << "Mean Detection efficiency:  " << meanEffTotal * 100 << "%" << std::endl;
    std::cout << "Cherenkov/Scintillation:    " << (meanScint > 0 ? meanCher/meanScint : 0) << std::endl;
    std::cout << "Total detected photons:     " << tPh->GetEntries() << std::endl;

    // Summary text on canvas
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.1, 0.9, "PbWO4 Optical Properties Summary");
    latex.DrawLatex(0.1, 0.82, Form("Total events: %d", nEvents));
    latex.DrawLatex(0.1, 0.74, Form("Mean Cherenkov photons/event: %.0f", meanCher));
    latex.DrawLatex(0.1, 0.66, Form("Mean Scintillation photons/event: %.0f", meanScint));
    latex.DrawLatex(0.1, 0.58, Form("Mean Detected photons/event: %.0f", meanDet));
    latex.DrawLatex(0.1, 0.50, Form("Mean Detection efficiency: %.1f%%", meanEffTotal * 100));
    latex.DrawLatex(0.1, 0.42, Form("Cherenkov/Scintillation ratio: %.2f",
                    meanScint > 0 ? meanCher/meanScint : 0));
    latex.DrawLatex(0.1, 0.34, Form("Total detected photons: %lld", tPh->GetEntries()));

    c4->SaveAs("plots_jpg/optical_summary.jpg");
    c4->Print("optical_results.pdf");

    // Close multi-page PDF
    cPdf->Print("optical_results.pdf]");

    file->Close();

    std::cout << "\n=== Output files generated ===" << std::endl;
    std::cout << "  optical_results.pdf                  - Complete Multi-page PDF report" << std::endl;
    std::cout << "  plots_jpg/production/optical_production.jpg - JPG image" << std::endl;
    std::cout << "  plots_jpg/spectrum/optical_spectrum.jpg   - JPG image" << std::endl;
    std::cout << "  plots_jpg/detection/optical_detection.jpg - JPG image" << std::endl;
    std::cout << "  plots_jpg/optical_summary.jpg             - JPG image" << std::endl;
}
