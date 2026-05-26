void pion_dual_readout_analysis(const char* filename = "output/output_thread-2.root")
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
        
        if (lastEvent != -1 && fEvent < lastEvent - 100) {
            currentRun++;
        }
        lastEvent = fEvent;
        
        int uniqueEvent = currentRun * 1000000 + fEvent;
        sumE[uniqueEvent] += fEdep;
        sumS[uniqueEvent] += fScint;
        sumC[uniqueEvent] += fCher;
    }

    std::ofstream out("pion_event_data.txt");
    out << "Event E S C\n";

    gStyle->SetOptStat(1111);

    // Full range for pion showers (hadronic: much more energy deposited + wider fluctuations)
    TH1D *hE   = new TH1D("hE",   "Energy Deposit (#pi^{-});E_{dep} [MeV];Events", 100, 0, 22000);
    TH1D *hS   = new TH1D("hS",   "Scintillation Signal (#pi^{-});S [a.u.];Events", 100, 0, 2.2e8);
    TH1D *hC   = new TH1D("hC",   "Cherenkov Signal (#pi^{-});C [a.u.];Events", 100, 0, 80000);
    TH1D *hSum = new TH1D("hSum", "Combined Signal S+C (#pi^{-});S+C [a.u.];Events", 100, 0, 2.2e8);

    TH2D *hE_S   = new TH2D("hE_S",   "S vs E_{dep} (#pi^{-});E_{dep} [MeV];S [a.u.]", 100, 0, 22000, 100, 0, 2.2e8);
    TH2D *hE_C   = new TH2D("hE_C",   "C vs E_{dep} (#pi^{-});E_{dep} [MeV];C [a.u.]", 100, 0, 22000, 100, 0, 80000);
    TH2D *hSC    = new TH2D("hSC",    "C vs S (#pi^{-});S [a.u.];C [a.u.]", 100, 0, 2.2e8, 100, 0, 80000);
    TH2D *hE_sum = new TH2D("hE_sum", "S+C vs E_{dep} (#pi^{-});E_{dep} [MeV];S+C [a.u.]", 100, 0, 22000, 100, 0, 2.2e8);

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

    TCanvas *c1 = new TCanvas("c1","Pion Distributions",1200,800);
    c1->Divide(2,2);
    c1->cd(1); hE->Draw();
    c1->cd(2); hS->Draw();
    c1->cd(3); hC->Draw();
    c1->cd(4); hSum->Draw();

    TCanvas *c2 = new TCanvas("c2","Pion Correlations",1200,800);
    c2->Divide(2,2);
    c2->cd(1); hE_S->Draw("COLZ");
    c2->cd(2); hE_C->Draw("COLZ");
    c2->cd(3); hSC->Draw("COLZ");
    c2->cd(4); hE_sum->Draw("COLZ");

    c1->SaveAs("pion_distributions.png");
    c2->SaveAs("pion_correlations.png");
}
