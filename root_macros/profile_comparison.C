/*
 * profile_comparison.C
 * Compares longitudinal and lateral shower profiles for pions vs muons
 * in a HOMOGENEOUS PbWO4 calorimeter (7x7x24, 10mm per layer).
 */
void profile_comparison(
    const char* pionFile = "pion_output/pion_homo.root",
    const char* muonFile = "muon_output/muon_homo.root")
{
    auto decodeID = [](int detID, int &layer, int &col, int &row) {
        layer = (detID / 10000) - 100;
        int local = detID % 10000;
        col = local / 100;
        row = local % 100;
    };

    const int nLayers = 24;
    const int nCols = 7;
    const int nRows = 7;
    const double layerThickness = 1.0; // cm (10 mm PbWO4)
    const double cellSize = 2.2; // cm

    // ============================================================
    // PROCESS PION FILE
    // ============================================================
    TFile *fPi = TFile::Open(pionFile);
    TTree *tPi = (TTree*)fPi->Get("Hits");

    Int_t fEvent; Double_t fEdep; Int_t fDetectorID;
    tPi->SetBranchAddress("fEvent", &fEvent);
    tPi->SetBranchAddress("fEdep", &fEdep);
    tPi->SetBranchAddress("fDetectorID", &fDetectorID);

    double piLongProfile[24] = {};
    double piLatProfile[7][7] = {};
    int piNevents = 0;
    std::map<int, double> piEventEnergy;
    std::map<int, std::map<int, double>> piRunEventEnergy;

    int lastEvt = -1, currentRun = 0;
    for (Long64_t i = 0; i < tPi->GetEntries(); i++) {
        tPi->GetEntry(i);
        if (lastEvt != -1 && fEvent < lastEvt - 100) currentRun++;
        lastEvt = fEvent;

        int layer, col, row;
        decodeID(fDetectorID, layer, col, row);
        if (layer < 0 || layer >= nLayers) continue;

        int uniqueEvt = currentRun * 1000000 + fEvent;
        piRunEventEnergy[currentRun][uniqueEvt] += fEdep;

        if (currentRun == 3) { // 10 GeV run
            piLongProfile[layer] += fEdep;
            if (col >= 0 && col < nCols && row >= 0 && row < nRows)
                piLatProfile[col][row] += fEdep;
            piEventEnergy[uniqueEvt] += fEdep;
        }
    }
    piNevents = piEventEnergy.size();
    if (piNevents > 0) {
        for (int l = 0; l < nLayers; l++) piLongProfile[l] /= piNevents;
        for (int c = 0; c < nCols; c++)
            for (int r = 0; r < nRows; r++)
                piLatProfile[c][r] /= piNevents;
    }

    // ============================================================
    // PROCESS MUON FILE
    // ============================================================
    TFile *fMu = TFile::Open(muonFile);
    TTree *tMu = (TTree*)fMu->Get("Hits");
    tMu->SetBranchAddress("fEvent", &fEvent);
    tMu->SetBranchAddress("fEdep", &fEdep);
    tMu->SetBranchAddress("fDetectorID", &fDetectorID);

    double muLongProfile[24] = {};
    double muLatProfile[7][7] = {};
    int muNevents = 0;
    std::map<int, double> muEventEnergy;
    std::map<int, std::map<int, double>> muRunEventEnergy;

    lastEvt = -1; currentRun = 0;
    for (Long64_t i = 0; i < tMu->GetEntries(); i++) {
        tMu->GetEntry(i);
        if (lastEvt != -1 && fEvent < lastEvt - 100) currentRun++;
        lastEvt = fEvent;

        int layer, col, row;
        decodeID(fDetectorID, layer, col, row);
        if (layer < 0 || layer >= nLayers) continue;

        int uniqueEvt = currentRun * 1000000 + fEvent;
        muRunEventEnergy[currentRun][uniqueEvt] += fEdep;

        if (currentRun == 3) { // 10 GeV
            muLongProfile[layer] += fEdep;
            if (col >= 0 && col < nCols && row >= 0 && row < nRows)
                muLatProfile[col][row] += fEdep;
            muEventEnergy[uniqueEvt] += fEdep;
        }
    }
    muNevents = muEventEnergy.size();
    if (muNevents > 0) {
        for (int l = 0; l < nLayers; l++) muLongProfile[l] /= muNevents;
        for (int c = 0; c < nCols; c++)
            for (int r = 0; r < nRows; r++)
                muLatProfile[c][r] /= muNevents;
    }

    printf("Pion 10 GeV events: %d, Muon 10 GeV events: %d\n", piNevents, muNevents);

    // ============================================================
    // PLOT 1: Longitudinal Profile Comparison (10 GeV)
    // ============================================================
    gStyle->SetOptStat(0);

    // X-axis in radiation lengths (X0 for PbWO4 = 0.89 cm, layer = 1.0 cm)
    double X0_PbWO4 = 0.89; // cm
    double layerInX0 = layerThickness / X0_PbWO4; // ~1.12 X0 per layer

    TH1D *hPiLong = new TH1D("hPiLong", "", nLayers, 0, nLayers * layerInX0);
    TH1D *hMuLong = new TH1D("hMuLong", "", nLayers, 0, nLayers * layerInX0);
    for (int l = 0; l < nLayers; l++) {
        hPiLong->SetBinContent(l+1, piLongProfile[l]);
        hMuLong->SetBinContent(l+1, muLongProfile[l]);
    }

    TCanvas *c1 = new TCanvas("c1", "Longitudinal Profile", 900, 600);
    c1->SetGrid();
    hPiLong->SetTitle("Longitudinal Profile: 10 GeV in PbWO_{4};Depth [X_{0}];Mean E_{dep} per layer [MeV]");
    hPiLong->SetLineColor(kRed+1);
    hPiLong->SetLineWidth(2);
    hPiLong->SetFillColorAlpha(kRed, 0.15);
    hPiLong->SetFillStyle(3004);
    hMuLong->SetLineColor(kBlue+1);
    hMuLong->SetLineWidth(2);
    hMuLong->SetFillColorAlpha(kBlue, 0.15);
    hMuLong->SetFillStyle(3005);

    double maxY = std::max(hPiLong->GetMaximum(), hMuLong->GetMaximum());
    hPiLong->GetYaxis()->SetRangeUser(0, maxY * 1.3);
    hPiLong->Draw("HIST");
    hMuLong->Draw("HIST SAME");

    TLegend *leg1 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg1->AddEntry(hPiLong, "#pi^{-} 10 GeV", "lf");
    leg1->AddEntry(hMuLong, "#mu^{+} 10 GeV", "lf");
    leg1->Draw();
    c1->SaveAs("longitudinal_comparison.png");

    // ============================================================
    // PLOT 2: Lateral Profile (radial from center) — 10 GeV
    // ============================================================
    const int centerCol = 3, centerRow = 3;

    const int nRadBins = 6;
    double piRadial[6] = {};
    double muRadial[6] = {};

    for (int c = 0; c < nCols; c++) {
        for (int r = 0; r < nRows; r++) {
            double dx = (c - centerCol) * cellSize;
            double dy = (r - centerRow) * cellSize;
            double dist = sqrt(dx*dx + dy*dy);
            int bin = (int)(dist / cellSize);
            if (bin >= nRadBins) bin = nRadBins - 1;
            piRadial[bin] += piLatProfile[c][r];
            muRadial[bin] += muLatProfile[c][r];
        }
    }

    TH1D *hPiRad = new TH1D("hPiRad", "", nRadBins, 0, nRadBins * cellSize);
    TH1D *hMuRad = new TH1D("hMuRad", "", nRadBins, 0, nRadBins * cellSize);
    for (int b = 0; b < nRadBins; b++) {
        double rInner = b * cellSize;
        double rOuter = (b+1) * cellSize;
        double area = M_PI * (rOuter*rOuter - rInner*rInner);
        hPiRad->SetBinContent(b+1, piRadial[b] / (area > 0 ? area : 1));
        hMuRad->SetBinContent(b+1, muRadial[b] / (area > 0 ? area : 1));
    }

    TCanvas *c2 = new TCanvas("c2", "Lateral Profile", 900, 600);
    c2->SetGrid();
    c2->SetLogy();
    hPiRad->SetTitle("Lateral Profile: 10 GeV in PbWO_{4};Radial Distance [cm];Mean E_{dep} / Area [MeV/cm^{2}]");
    hPiRad->SetLineColor(kRed+1);  hPiRad->SetLineWidth(2);
    hPiRad->SetMarkerColor(kRed+1); hPiRad->SetMarkerStyle(20); hPiRad->SetMarkerSize(1.3);
    hMuRad->SetLineColor(kBlue+1); hMuRad->SetLineWidth(2);
    hMuRad->SetMarkerColor(kBlue+1); hMuRad->SetMarkerStyle(21); hMuRad->SetMarkerSize(1.3);

    hPiRad->Draw("HIST P");
    hMuRad->Draw("HIST P SAME");

    TLegend *leg2 = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg2->AddEntry(hPiRad, "#pi^{-} 10 GeV", "lp");
    leg2->AddEntry(hMuRad, "#mu^{+} 10 GeV", "lp");
    leg2->Draw();
    c2->SaveAs("lateral_comparison.png");

    // ============================================================
    // PLOT 3: Lateral heatmaps side by side
    // ============================================================
    TH2D *hPiHeat = new TH2D("hPiHeat", "#pi^{-} 10 GeV Lateral Map;Column;Row",
                              nCols, -0.5, nCols-0.5, nRows, -0.5, nRows-0.5);
    TH2D *hMuHeat = new TH2D("hMuHeat", "#mu^{+} 10 GeV Lateral Map;Column;Row",
                              nCols, -0.5, nCols-0.5, nRows, -0.5, nRows-0.5);
    for (int c = 0; c < nCols; c++) {
        for (int r = 0; r < nRows; r++) {
            hPiHeat->SetBinContent(c+1, r+1, piLatProfile[c][r]);
            hMuHeat->SetBinContent(c+1, r+1, muLatProfile[c][r]);
        }
    }

    TCanvas *c3 = new TCanvas("c3", "Lateral Heatmaps", 1200, 500);
    c3->Divide(2, 1);
    c3->cd(1); gPad->SetRightMargin(0.15); hPiHeat->Draw("COLZ TEXT");
    c3->cd(2); gPad->SetRightMargin(0.15); hMuHeat->Draw("COLZ TEXT");
    c3->SaveAs("lateral_heatmaps.png");

    // ============================================================
    // PLOT 4: Energy deposit distributions at 10 GeV
    // ============================================================
    TH1D *hPiEdep = new TH1D("hPiEdep", "Total E_{dep} at 10 GeV (PbWO_{4});E_{dep} [MeV];Events",
                              100, 0, 11000);
    TH1D *hMuEdep = new TH1D("hMuEdep", "", 100, 0, 11000);

    for (auto &kv : piEventEnergy) hPiEdep->Fill(kv.second);
    for (auto &kv : muEventEnergy) hMuEdep->Fill(kv.second);

    TCanvas *c4 = new TCanvas("c4", "Edep Distributions", 900, 600);
    c4->SetGrid();
    hPiEdep->SetLineColor(kRed+1); hPiEdep->SetLineWidth(2);
    hMuEdep->SetLineColor(kBlue+1); hMuEdep->SetLineWidth(2);

    double max4 = std::max(hPiEdep->GetMaximum(), hMuEdep->GetMaximum());
    hPiEdep->GetYaxis()->SetRangeUser(0, max4 * 1.3);
    hPiEdep->Draw("HIST");
    hMuEdep->Draw("HIST SAME");

    TLegend *leg4 = new TLegend(0.45, 0.7, 0.88, 0.88);
    leg4->AddEntry(hPiEdep, Form("#pi^{-} 10 GeV: #mu=%.0f, #sigma=%.0f MeV",
                   hPiEdep->GetMean(), hPiEdep->GetStdDev()), "l");
    leg4->AddEntry(hMuEdep, Form("#mu^{+} 10 GeV: #mu=%.1f, #sigma=%.1f MeV",
                   hMuEdep->GetMean(), hMuEdep->GetStdDev()), "l");
    leg4->Draw();
    c4->SaveAs("edep_comparison.png");

    // ============================================================
    // PLOT 5: Energy resolution vs beam energy
    // ============================================================
    double beamE[5] = {1, 2, 5, 10, 20};
    double piRes[5], muRes[5], piMeanArr[5], muMeanArr[5];

    for (int run = 0; run < 5; run++) {
        if (piRunEventEnergy.count(run)) {
            double sum = 0, sum2 = 0; int n = 0;
            for (auto &kv : piRunEventEnergy[run]) {
                sum += kv.second; sum2 += kv.second * kv.second; n++;
            }
            double mean = sum / n;
            double sigma = sqrt(sum2/n - mean*mean);
            piMeanArr[run] = mean;
            piRes[run] = sigma / mean;
        }
        if (muRunEventEnergy.count(run)) {
            double sum = 0, sum2 = 0; int n = 0;
            for (auto &kv : muRunEventEnergy[run]) {
                sum += kv.second; sum2 += kv.second * kv.second; n++;
            }
            double mean = sum / n;
            double sigma = sqrt(sum2/n - mean*mean);
            muMeanArr[run] = mean;
            muRes[run] = sigma / mean;
        }
    }

    // Find max res for axis
    double maxRes = 0;
    for (int i = 0; i < 5; i++) {
        if (piRes[i] > maxRes) maxRes = piRes[i];
        if (muRes[i] > maxRes) maxRes = muRes[i];
    }

    TGraph *gPiRes = new TGraph(5, beamE, piRes);
    TGraph *gMuRes = new TGraph(5, beamE, muRes);

    TCanvas *c5 = new TCanvas("c5", "Energy Resolution", 900, 600);
    c5->SetGrid();
    c5->SetLogx();

    gPiRes->SetTitle("Energy Resolution in PbWO_{4}: #sigma_{E}/E vs Beam Energy;E_{beam} [GeV];#sigma_{E} / E");
    gPiRes->SetMarkerStyle(20); gPiRes->SetMarkerColor(kRed+1); gPiRes->SetLineColor(kRed+1);
    gPiRes->SetMarkerSize(1.5); gPiRes->SetLineWidth(2);
    gMuRes->SetMarkerStyle(21); gMuRes->SetMarkerColor(kBlue+1); gMuRes->SetLineColor(kBlue+1);
    gMuRes->SetMarkerSize(1.5); gMuRes->SetLineWidth(2);

    gPiRes->GetYaxis()->SetRangeUser(0, maxRes * 1.2);
    gPiRes->GetXaxis()->SetLimits(0.8, 25);
    gPiRes->Draw("APL");
    gMuRes->Draw("PL SAME");

    TLatex latex;
    latex.SetTextSize(0.03);
    for (int i = 0; i < 5; i++) {
        latex.SetTextColor(kRed+1);
        latex.DrawLatex(beamE[i]*1.08, piRes[i]+0.005, Form("%.1f%%", piRes[i]*100));
        latex.SetTextColor(kBlue+1);
        latex.DrawLatex(beamE[i]*1.08, muRes[i]+0.005, Form("%.1f%%", muRes[i]*100));
    }

    TLegend *leg5 = new TLegend(0.15, 0.75, 0.50, 0.88);
    leg5->AddEntry(gPiRes, "#pi^{-} (hadronic)", "lp");
    leg5->AddEntry(gMuRes, "#mu^{+} (MIP)", "lp");
    leg5->Draw();
    c5->SaveAs("resolution_comparison.png");

    // ============================================================
    printf("\n=== SUMMARY TABLE (Homogeneous PbWO4) ===\n");
    printf("%-7s | %-12s %-12s %-12s | %-12s %-12s %-12s\n",
           "E_beam", "pi <Edep>", "pi sigma", "pi sig/E",
           "mu <Edep>", "mu sigma", "mu sig/E");
    printf("--------|----------------------------------------|----------------------------------------\n");
    for (int r = 0; r < 5; r++) {
        double piSig = piRes[r] * piMeanArr[r];
        double muSig = muRes[r] * muMeanArr[r];
        printf("%-4.0f GeV| %8.1f MeV %8.1f MeV %10.4f   | %8.1f MeV %8.1f MeV %10.4f\n",
               beamE[r], piMeanArr[r], piSig, piRes[r], muMeanArr[r], muSig, muRes[r]);
    }
}
