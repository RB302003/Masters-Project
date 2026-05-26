void plotLongitudinalWithPeaks_fromDetectorID() {
    gStyle->SetOptStat(0);

    const char* fileName = "/home/rudradeb/sim/output/output_thread-1.root";
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "❌ Could not open file: " << fileName << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Hits");
    if (!tree) {
        std::cerr << "❌ Tree 'Hits' not found in file: " << fileName << std::endl;
        file->Close();
        return;
    }

    Int_t fDetectorID;
    Double_t fEdep;
    tree->SetBranchAddress("fDetectorID", &fDetectorID);
    tree->SetBranchAddress("fEdep", &fEdep);

    const int nLayers = 24;
    double sumE[nLayers] = {0};
    double maxE[nLayers] = {0};
    int countE[nLayers] = {0};

    Long64_t nentries = tree->GetEntries();
    std::cout << "Total entries: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        int layer = (fDetectorID / 10000) - 100;  // active detIDs have (layer+100)*10000

        // If your layer numbering starts from 1 (e.g. 1..24) and you want 0-based,
        // uncomment the next line:
        // layer = layer - 1;

        if (layer >= 0 && layer < nLayers) {
            sumE[layer] += fEdep;
            countE[layer]++;
            if (fEdep > maxE[layer]) maxE[layer] = fEdep;
        }
    }
    file->Close();

    // Debugging output
    std::cout << "Entries per layer:" << std::endl;
    for (int l = 0; l < nLayers; ++l) {
        std::cout << "Layer " << l << ": " << countE[l] << " entries" << std::endl;
    }

    // Create graphs
    TCanvas *c = new TCanvas("c", "Longitudinal Shower Profiles", 1200, 600);
    c->Divide(2,1);

    // Total Energy per Layer
    c->cd(1);
    TGraph *grSum = new TGraph(nLayers);
    for (int l = 0; l < nLayers; ++l)
        grSum->SetPoint(l, l+1, sumE[l]);
    grSum->SetTitle("Longitudinal Profile - Total Energy");
    grSum->GetXaxis()->SetTitle("Layer");
    grSum->GetYaxis()->SetTitle("Deposited Energy (MeV)");
    grSum->SetMarkerStyle(20);
    grSum->SetMarkerColor(kRed);
    grSum->SetLineColor(kRed);
    grSum->Draw("APL");

    // Peak Energy per Layer
    c->cd(2);
    TGraph *grMax = new TGraph(nLayers);
    for (int l = 0; l < nLayers; ++l)
        grMax->SetPoint(l, l+1, maxE[l]);
    grMax->SetTitle("Longitudinal Profile - Peak Energy per Layer");
    grMax->GetXaxis()->SetTitle("Layer");
    grMax->GetYaxis()->SetTitle("Max Deposited Energy (MeV)");
    grMax->SetMarkerStyle(21);
    grMax->SetMarkerColor(kBlue);
    grMax->SetLineColor(kBlue);
    grMax->Draw("APL");

    c->Update();
}

