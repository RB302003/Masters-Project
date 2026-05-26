void dual_readout_analysis(const char* filename = "CherenkovOutput/output_thread-2.root")
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }

    TTree *tree = (TTree*)f->Get("Hits");
    if (!tree) { printf("No 'Hits' tree found\n"); return; }

    Int_t fEvent;
    Double_t fEdep = 0;
    Double_t fScint = 0;
    Double_t fCher = 0;

    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("fEdep", &fEdep);
    tree->SetBranchAddress("fScint", &fScint);
    tree->SetBranchAddress("fCher", &fCher);

    std::map<int,double> sumE, sumS, sumC;

    Long64_t nEntries = tree->GetEntries();

    int currentRun = 0;
    int lastEvent = -1;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Detect when fEvent resets (new /run/beamOn)
        if (lastEvent != -1 && fEvent < lastEvent - 100) {
            currentRun++;
        }
        lastEvent = fEvent;
        
        int uniqueEvent = currentRun * 1000000 + fEvent;
        sumE[uniqueEvent] += fEdep;
        sumS[uniqueEvent] += fScint;
        sumC[uniqueEvent] += fCher;
    }

    std::ofstream out("event_data.txt");
    out << "Event E S C\n";

    gSystem->Exec("mkdir -p plots_jpg/distributions");
    gSystem->Exec("mkdir -p plots_jpg/correlations");

    TCanvas* cPdf = new TCanvas("cPdf", "PDF Wrapper", 1200, 800);
    cPdf->Print("dual_readout_results.pdf[");

    gStyle->SetOptStat(1111);

    TH1D *hE   = new TH1D("hE",   "Energy Deposit;E_{dep} [MeV];Events", 100, 0, 1200);
    TH1D *hS   = new TH1D("hS",   "Scintillation Signal;S [a.u.];Events", 100, 0, 1.2e7);
    TH1D *hC   = new TH1D("hC",   "Cherenkov Signal;C [a.u.];Events", 100, 0, 12000);
    TH1D *hSum = new TH1D("hSum", "Combined Signal (S+C);S+C [a.u.];Events", 100, 0, 1.2e7);

    TH2D *hE_S   = new TH2D("hE_S",   "S vs E_{dep};E_{dep} [MeV];S [a.u.]", 100, 0, 1200, 100, 0, 1.2e7);
    TH2D *hE_C   = new TH2D("hE_C",   "C vs E_{dep};E_{dep} [MeV];C [a.u.]", 100, 0, 1200, 100, 0, 12000);
    TH2D *hSC    = new TH2D("hSC",    "C vs S;S [a.u.];C [a.u.]", 100, 0, 1.2e7, 100, 0, 12000);
    TH2D *hE_sum = new TH2D("hE_sum", "S+C vs E_{dep};E_{dep} [MeV];S+C [a.u.]", 100, 0, 1200, 100, 0, 1.2e7);

    for (auto &kv : sumE) {
        int evt = kv.first;

        double E = sumE[evt];
        double S = sumS[evt];
        double C = sumC[evt];

        out << evt << " " << E << " " << S << " " << C << "\n";

        hE->Fill(E);
        hS->Fill(S);
        hC->Fill(C);
        hSum->Fill(S + C);

        hE_S->Fill(E,S);
        hE_C->Fill(E,C);
        hSC->Fill(S,C);
        hE_sum->Fill(E,S+C);
    }

    out.close();

    TCanvas *c1 = new TCanvas("c1","Distributions",1200,800);
    c1->Divide(2,2);
    c1->cd(1); hE->Draw();
    c1->cd(2); hS->Draw();
    c1->cd(3); hC->Draw();
    c1->cd(4); hSum->Draw();

    c1->SaveAs("plots_jpg/distributions/distributions.jpg");
    c1->Print("dual_readout_results.pdf");

    TCanvas *c2 = new TCanvas("c2","Correlations",1200,800);
    c2->Divide(2,2);
    c2->cd(1); hE_S->Draw("COLZ");
    c2->cd(2); hE_C->Draw("COLZ");
    c2->cd(3); hSC->Draw("COLZ");
    c2->cd(4); hE_sum->Draw("COLZ");

    c2->SaveAs("plots_jpg/correlations/correlations.jpg");
    c2->Print("dual_readout_results.pdf");

    cPdf->Print("dual_readout_results.pdf]");
    
    std::cout << "\n=== Output files generated ===" << std::endl;
    std::cout << "  dual_readout_results.pdf                      - Multi-page PDF report" << std::endl;
    std::cout << "  plots_jpg/distributions/distributions.jpg     - JPG image" << std::endl;
    std::cout << "  plots_jpg/correlations/correlations.jpg       - JPG image (Exactly matching your requested diagram)" << std::endl;
}
