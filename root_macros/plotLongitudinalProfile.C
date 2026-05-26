// ============================================================================
//  Longitudinal Shower Profile Plotter
//  -----------------------------------
//  This program reads energy deposition data from your ROOT simulation files
//  and plots how the shower develops *along the depth* of the calorimeter.
//
//  ✅ X-axis = Calorimeter layers (1 to 24)
//  ✅ Y-axis = Energy deposited per layer
//  ✅ One graph each for 10, 20, 50, 100 GeV runs
//  ✅ Uses output files from your /home/rudradeb/sim/finaloutput/
//
//  Author: <Your name or leave blank>
// ============================================================================

void plotLongitudinalProfile(Bool_t normalizeToPeak = kTRUE, Bool_t perEventAverage = kFALSE)
{
    // ----------------------------- Settings ---------------------------------
    // normalizeToPeak = kTRUE  → Each graph is normalized by its highest layer value.
    //                             This lets you compare *shape only* between energies.
    //
    // perEventAverage = kFALSE → If TRUE, divides energy by total number of events
    //                             to find average energy per layer per event.
    // -------------------------------------------------------------------------

    gStyle->SetOptStat(0); // Turn off default statistics box from ROOT

    // ----------------------------------------------------------------------------
    //  File mapping: associate each run energy with its ROOT file
    //  You said:
    //      10 GeV  → output_thread-1.root
    //      20 GeV  → output_thread-2.root
    //      50 GeV  → output_thread-3.root
    //      100 GeV → output_thread-4.root
    // ----------------------------------------------------------------------------
    std::map<int, std::string> files = {
        {1,  "/home/rudradeb/sim/muonoutput/output_thread-1.root"},
        {2,  "/home/rudradeb/sim/muonoutput/output_thread-2.root"},
        {5,  "/home/rudradeb/sim/muonoutput/output_thread-3.root"},
        {10, "/home/rudradeb/sim/muonoutput/output_thread-4.root"}
    };

    // ----------------------------------------------------------------------------
    // Helper to extract the *layer number* from fDetectorID
    // Your ID format = layer*10000 + col*100 + row
    // ----------------------------------------------------------------------------
    auto decodeLayer = [](int detID) {
        return (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    };

    // The detector has 24 layers
    const int nLayers = 24;

    // Colors for plotting 4 runs
    const int colors[4] = {kRed, kBlue, kGreen + 2, kMagenta};
    const int energies[4] = {1, 2, 5, 10};

    // -------------------- Create a 2×2 plotting canvas -----------------------
    TCanvas* c = new TCanvas("c", "Longitudinal Shower Profiles", 1200, 900);
    c->Divide(2, 2); // 4 graphs = 4 runs

    // --------------------------- Loop over each energy file ------------------
    for (int i = 0; i < 4; ++i)
    {
        int E = energies[i];
        TFile* file = TFile::Open(files[E].c_str());
        if (!file || file->IsZombie()) {
            printf("Error: Cannot open file for %d GeV run.\n", E);
            continue;
        }

        // Get the "Hits" TTree from ROOT file
        TTree* tree = (TTree*)file->Get("Hits");
        if (!tree) {
            printf("Error: 'Hits' tree not found in %s\n", files[E].c_str());
            continue;
        }

        // Variables connected to ROOT tree branches
        Int_t detID;
        Double_t edep;
        Int_t eventID;

        // Link them
        tree->SetBranchAddress("fDetectorID", &detID);
        tree->SetBranchAddress("fEdep", &edep);
        tree->SetBranchAddress("fEvent", &eventID);

        // We will collect energy in each layer here
        std::vector<double> layerEnergy(nLayers, 0.0);

        // If user wants average per event, count unique events
        std::set<int> eventCount;

        // ---------------- Loop over all hits in the file --------------------
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t j = 0; j < nEntries; j++) {
            tree->GetEntry(j);

            // Ignore non-contributing hits
            if (edep <= 0) continue;

            // Extract layer index from Detector ID
            int layer = decodeLayer(detID);
            if (layer >= 0 && layer < nLayers) {  // valid 0-based layer index
                // Add energy to the corresponding layer
                layerEnergy[layer] += edep;
                eventCount.insert(eventID);
            }
        }

        // -------- Optional: convert total to per-event average --------------
        if (perEventAverage && !eventCount.empty()) {
            for (double &E : layerEnergy) E /= eventCount.size();
        }

        // -------- Optional: Normalize each profile to its own maximum --------
        if (normalizeToPeak) {
            double maxVal = *max_element(layerEnergy.begin(), layerEnergy.end());
            if (maxVal > 0)
                for (double &E : layerEnergy) E /= maxVal;
        }

        // ------------------------ Plot the results ---------------------------
        c->cd(i + 1); // Go to the next panel (total 4)

        TGraph* gr = new TGraph(nLayers);
        for (int L = 0; L < nLayers; ++L)
            gr->SetPoint(L, L + 1, layerEnergy[L]); // (Layer, Energy)

        gr->SetTitle(Form("Longitudinal Profile - %d GeV", E));
        gr->GetXaxis()->SetTitle("Layer Number");
        gr->GetYaxis()->SetTitle(normalizeToPeak ?
            "Normalized Energy Deposition" :
            "Total Energy Deposition (MeV)");
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]);
        gr->Draw("ALP");

        file->Close();
    }

    c->Update();
}

