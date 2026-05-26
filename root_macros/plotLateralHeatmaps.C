// plotLateralHeatmaps_safe.C
// Safer variant of plotLateral_showerDev_realistic_fix.C
// Paste into ROOT and run: plotLateralHeatmaps_safe();

#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom3.h"

// decode detector ID format: L*10000 + col*100 + row
static inline void decodeDet(int detID, int &L, int &col, int &row) {
    L  = (detID / 10000) - 100;  // active detIDs have (layer+100)*10000
    int tmp = detID % 10000;
    col = tmp / 100;
    row = tmp % 100;
}

// small gaussian kernel builder
static void makeGaussKernel(double sigma, std::vector<std::vector<double>> &kern, int &R) {
    R = std::max(1, (int)ceil(3.0 * sigma));
    kern.assign(2*R+1, std::vector<double>(2*R+1,0.0));
    double sum=0;
    for(int i=-R;i<=R;i++) for(int j=-R;j<=R;j++){
        double v = std::exp(-0.5*(i*i + j*j)/(sigma*sigma));
        kern[i+R][j+R]=v; sum+=v;
    }
    for(int i=0;i<2*R+1;i++) for(int j=0;j<2*R+1;j++) kern[i][j]/=sum;
}

// bilinear sample from 2D grid 'g' of size Nx x Ny at continuous coords (u,v).
// returns 0 if (u,v) outside [0,Nx-1]x[0,Ny-1]
static double bilinearSample(const std::vector<std::vector<double>> &g, double u, double v) {
    int Nx = (int)g.size();
    if (Nx == 0) return 0.0;
    int Ny = (int)g[0].size();
    if (u < 0.0 || v < 0.0 || u > (double)(Nx-1) || v > (double)(Ny-1)) return 0.0;
    int i = (int)floor(u);
    int j = (int)floor(v);
    double du = u - i;
    double dv = v - j;
    int i1 = std::min(i+1, Nx-1);
    int j1 = std::min(j+1, Ny-1);
    double v00 = g[i][j];
    double v10 = g[i1][j];
    double v01 = g[i][j1];
    double v11 = g[i1][j1];
    return (1-du)*(1-dv)*v00 + du*(1-dv)*v10 + (1-du)*dv*v01 + du*dv*v11;
}

void plotLateralHeatmaps(Bool_t perEventNormalize=false,
                         Double_t energyScale = 1e-3,   // default: MeV -> GeV
                         Double_t smoothSigma = 0.6,    // final mild smoothing
                         Int_t oversample = 4)          // subpixel oversampling factor
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(140);

    // --- files/energies (edit paths if needed) ---
    std::map<int,std::string> fileOf = {
        {1,  "/home/rudradeb/sim/muonoutput/output_thread-1.root"},
        {2,  "/home/rudradeb/sim/muonoutput/output_thread-2.root"},
        {5,  "/home/rudradeb/sim/muonoutput/output_thread-3.root"},
        {10, "/home/rudradeb/sim/muonoutput/output_thread-4.root"}
    };
    std::vector<int> energies = {1,2,5,10};
    // ----------------------------------------------

    const int nLayers = 24;
    const int outSide = 7;        // final heatmap size (7x7)
    const int eventSide = 11;     // working grid (7 centered in 11)
    const int eventCenter = (eventSide-1)/2;
    const int halfOut = (outSide-1)/2;

    // final convolution kernel
    std::vector<std::vector<double>> kernel; int kR=1; makeGaussKernel(smoothSigma, kernel, kR);

    for (int E : energies) {
        std::printf("\nProcessing %d GeV (scale=%.3g, oversample=%d)...\n", E, energyScale, oversample);
        auto it = fileOf.find(E);
        if (it == fileOf.end()) { Warning("plotLateralHeatmaps_safe","No filename mapped for %d GeV", E); continue; }
        std::string fname = it->second;

        TFile fin(fname.c_str(), "READ");
        if (!fin.IsOpen() || fin.IsZombie()) {
            Warning("plotLateralHeatmaps_safe","Cannot open %s", fname.c_str());
            if (fin.IsOpen()) fin.Close();
            continue;
        }

        TTree *tree = dynamic_cast<TTree*>(fin.Get("Hits"));
        if (!tree) {
            Warning("plotLateralHeatmaps_safe","No Hits TTree in %s", fname.c_str());
            fin.Close();
            continue;
        }

        Int_t detID=0, trID=0, evID=0; Double_t edep=0;
        tree->SetBranchAddress("fDetectorID",&detID);
        tree->SetBranchAddress("fEdep",&edep);
        // guard optional branches: only set if they exist
        if (tree->GetBranch("fTrackID")) tree->SetBranchAddress("fTrackID",&trID);
        if (tree->GetBranch("fEvent"))   tree->SetBranchAddress("fEvent",&evID);

        // 1) accumulate per-event per-layer 7x7 maps (column,row)
        std::map<Int_t, std::vector<std::vector<std::vector<double>>>> perEventLayer;
        Long64_t N = tree->GetEntries();
        for (Long64_t i=0;i<N;++i) {
            tree->GetEntry(i);
            if (edep <= 0) continue;
            int L, col, row; decodeDet(detID,L,col,row);
            if (L < 0 || L >= nLayers) continue;
            if (col < 0 || col >= outSide || row < 0 || row >= outSide) continue;
            if (perEventLayer.find(evID) == perEventLayer.end()) {
                perEventLayer[evID] = std::vector<std::vector<std::vector<double>>>(nLayers,
                    std::vector<std::vector<double>>(outSide, std::vector<double>(outSide, 0.0)));
            }
            // accumulate
            perEventLayer[evID][L][col][row] += edep * energyScale; // convert to GeV
        }
        fin.Close();

        if (perEventLayer.empty()) {
            Warning("plotLateralHeatmaps_safe","No events found for %d GeV (empty perEventLayer).", E);
            continue;
        }

        // 2) build Event objects (centroid computed on near-peak high-res map)
        struct Event {
            int peakL;
            double cx; // centroid x (column) in eventSide coords [0..eventSide-1], fractional
            double cy; // centroid y (row)
            std::vector<std::vector<std::vector<double>>> layers; // [nLayers][7][7]
        };
        std::vector<Event> events; events.reserve(perEventLayer.size());

        for (auto &pp : perEventLayer) {
            auto &layerCell = pp.second;
            // compute layer sums to find peak
            std::vector<double> sums(nLayers,0.0);
            for (int L=0; L<nLayers; ++L)
                for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) sums[L] += layerCell[L][c][r];
            // find peak layer safely
            auto itmax = std::max_element(sums.begin(), sums.end());
            if (itmax == sums.end()) continue;
            int peakL = (int)std::distance(sums.begin(), itmax);

            // near-peak band for centroiding
            int p0 = std::max(0, peakL-1);
            int p1 = std::min(nLayers-1, peakL+1);

            // build high-res map: eventSide * oversample
            int highSide = eventSide * oversample;
            if (highSide <= 0) continue;
            std::vector<std::vector<double>> highMap(highSide, std::vector<double>(highSide, 0.0));
            int offset = eventCenter - halfOut;
            int offsetHigh = offset * oversample;

            // deposit near-peak energy into high-res grid (spread each coarse cell into an oversample block)
            for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
                double val = 0;
                for (int L=p0; L<=p1; ++L) val += layerCell[L][c][r];
                if (val == 0.0) continue;
                for (int ix=0; ix<oversample; ++ix) for (int iy=0; iy<oversample; ++iy) {
                    int hx = offsetHigh + c*oversample + ix;
                    int hy = offsetHigh + r*oversample + iy;
                    if (hx>=0 && hx<highSide && hy>=0 && hy<highSide) highMap[hx][hy] += val / (oversample*oversample);
                }
            }

            // small 3x3 box blur to regularize centroid (cheap)
            std::vector<std::vector<double>> tmp(highSide, std::vector<double>(highSide,0.0));
            for (int i=0;i<highSide;++i) for (int j=0;j<highSide;++j) {
                double s=0; int cnt=0;
                for (int di=-1; di<=1; ++di) for (int dj=-1; dj<=1; ++dj) {
                    int ii=i+di, jj=j+dj;
                    if (ii<0||ii>=highSide||jj<0||jj>=highSide) continue;
                    s += highMap[ii][jj]; cnt++;
                }
                tmp[i][j] = (cnt>0) ? s/cnt : 0.0;
            }
            highMap.swap(tmp);

            // centroid in high-res coords
            double stot=0, sx=0, sy=0;
            for (int i=0;i<highSide;++i) for (int j=0;j<highSide;++j) {
                double v = highMap[i][j];
                stot += v;
                sx += v * (double)i;
                sy += v * (double)j;
            }
            if (stot <= 0) continue; // skip empty events
            double cx_high = sx / stot;
            double cy_high = sy / stot;
            // convert back to eventSide fractional coords
            double cx = cx_high / (double)oversample;
            double cy = cy_high / (double)oversample;

            Event ev;
            ev.peakL = peakL;
            ev.cx = cx;
            ev.cy = cy;
            ev.layers = layerCell; // copy
            events.push_back(std::move(ev));
        }
        perEventLayer.clear();

        if (events.empty()) {
            Warning("plotLateralHeatmaps_safe","No usable events after centroiding for %d GeV.", E);
            continue;
        }
        std::printf(" events used: %zu\n", events.size());

        // 3) accumulators for EARLY / NEAR-PEAK / LATE
        std::vector<std::vector<double>> accumE(outSide, std::vector<double>(outSide,0.0));
        std::vector<std::vector<double>> accumP(outSide, std::vector<double>(outSide,0.0));
        std::vector<std::vector<double>> accumL(outSide, std::vector<double>(outSide,0.0));
        int usedEvents = 0;

        for (auto &ev : events) {
            // per-cell bands (coarse 7x7)
            std::vector<std::vector<double>> bandE(outSide, std::vector<double>(outSide,0.0));
            std::vector<std::vector<double>> bandP(outSide, std::vector<double>(outSide,0.0));
            std::vector<std::vector<double>> bandL(outSide, std::vector<double>(outSide,0.0));

            int peak = ev.peakL;
            int p0 = std::max(0, peak-1);
            int p1 = std::min(nLayers-1, peak+1);

            for (int L=0; L<nLayers; ++L) {
                for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
                    double v = ev.layers[L][c][r];
                    if (L <= 4) bandE[c][r] += v;
                    else if (L >= p0 && L <= p1) bandP[c][r] += v;
                    else bandL[c][r] += v; // everything else to late (tail)
                }
            }

            // event total (coarse units)
            double Etot=0;
            for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) Etot += bandE[c][r] + bandP[c][r] + bandL[c][r];
            if (Etot <= 0) continue;

            // create high-res band maps by spreading each coarse cell into oversample blocks
            int highSide = eventSide * oversample;
            int offset = eventCenter - halfOut;
            int offsetHigh = offset * oversample;
            std::vector<std::vector<double>> hE(highSide, std::vector<double>(highSide,0.0));
            std::vector<std::vector<double>> hP(highSide, std::vector<double>(highSide,0.0));
            std::vector<std::vector<double>> hL(highSide, std::vector<double>(highSide,0.0));

            for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
                double vE = bandE[c][r];
                double vP = bandP[c][r];
                double vL = bandL[c][r];
                if (vE==0 && vP==0 && vL==0) continue;
                for (int ix=0; ix<oversample; ++ix) for (int iy=0; iy<oversample; ++iy) {
                    int hx = offsetHigh + c*oversample + ix;
                    int hy = offsetHigh + r*oversample + iy;
                    if (hx>=0 && hx<highSide && hy>=0 && hy<highSide) {
                        hE[hx][hy] += vE / (oversample*oversample);
                        hP[hx][hy] += vP / (oversample*oversample);
                        hL[hx][hy] += vL / (oversample*oversample);
                    }
                }
            }

            // compute shift so that event centroid lands at eventCenter
            double cx_high = ev.cx * (double)oversample;
            double cy_high = ev.cy * (double)oversample;
            double target_high = eventCenter * (double)oversample;
            double shift_x = cx_high - target_high;
            double shift_y = cy_high - target_high;

            // sample back to coarse 7x7 by sampling high-res maps at shifted coords
            for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
                // coarse cell center in eventSide coords
                double tx_coarse = (double)(c + offset);
                double ty_coarse = (double)(r + offset);
                // high-res center of that coarse block
                double tx_high = tx_coarse * (double)oversample + 0.5*(oversample-1);
                double ty_high = ty_coarse * (double)oversample + 0.5*(oversample-1);
                // source sample point = target - shift (we move event to center)
                double srcx = tx_high - shift_x;
                double srcy = ty_high - shift_y;

                double valE = bilinearSample(hE, srcx, srcy);
                double valP = bilinearSample(hP, srcx, srcy);
                double valL = bilinearSample(hL, srcx, srcy);

                double outE = perEventNormalize ? (valE / Etot) : valE;
                double outP = perEventNormalize ? (valP / Etot) : valP;
                double outL = perEventNormalize ? (valL / Etot) : valL;

                accumE[c][r] += outE;
                accumP[c][r] += outP;
                accumL[c][r] += outL;
            }

            usedEvents++;
        } // events

        if (usedEvents == 0) { Warning("plotLateralHeatmaps_safe","No used events for %d", E); continue; }
        // average
        for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
            accumE[c][r] /= usedEvents;
            accumP[c][r] /= usedEvents;
            accumL[c][r] /= usedEvents;
        }

        // final mild smoothing (convolution)
        auto convolve = [&](const std::vector<std::vector<double>> &in)->std::vector<std::vector<double>> {
            std::vector<std::vector<double>> out(outSide, std::vector<double>(outSide,0.0));
            for (int x=0;x<outSide;++x) for (int y=0;y<outSide;++y) {
                double s=0;
                for (int i=-kR;i<=kR;i++) for (int j=-kR;j<=kR;j++) {
                    int xx = x + i; int yy = y + j;
                    if (xx<0||xx>=outSide||yy<0||yy>=outSide) continue;
                    s += in[xx][yy] * kernel[i+kR][j+kR];
                }
                out[x][y] = s;
            }
            return out;
        };

        auto eE = convolve(accumE);
        auto eP = convolve(accumP);
        auto eL = convolve(accumL);

        // fill TH2 (x=column, y=row)
        std::string nameE = Form("hEARLY_%d",E);
        std::string nameP = Form("hPEAK_%d",E);
        std::string nameL = Form("hLATE_%d",E);
        TH2D *hE = new TH2D(nameE.c_str(),"EARLY (layers 0-4)", outSide, -0.5, outSide-0.5, outSide, -0.5, outSide-0.5);
        TH2D *hP = new TH2D(nameP.c_str(),"NEAR-PEAK (peak±1)", outSide, -0.5, outSide-0.5, outSide, -0.5, outSide-0.5);
        TH2D *hL = new TH2D(nameL.c_str(),"LATE (tail)", outSide, -0.5, outSide-0.5, outSide, -0.5, outSide-0.5);

        for (int c=0;c<outSide;++c) for (int r=0;r<outSide;++r) {
            hE->SetBinContent(c+1, r+1, eE[c][r]);
            hP->SetBinContent(c+1, r+1, eP[c][r]);
            hL->SetBinContent(c+1, r+1, eL[c][r]);
        }

        const char *zTitle = perEventNormalize ? "Avg fraction / event" : "Energy (GeV)";
        hE->GetZaxis()->SetTitle(zTitle);
        hP->GetZaxis()->SetTitle(zTitle);
        hL->GetZaxis()->SetTitle(zTitle);
        hE->GetXaxis()->SetTitle("Column"); hE->GetYaxis()->SetTitle("Row");
        hP->GetXaxis()->SetTitle("Column"); hP->GetYaxis()->SetTitle("Row");
        hL->GetXaxis()->SetTitle("Column"); hL->GetYaxis()->SetTitle("Row");

        // draw
        TCanvas *c = new TCanvas(Form("latdev_%d",E), Form("Lateral development %d GeV", E), 1260, 420);
        c->Divide(3,1);
        gStyle->SetPaintTextFormat(".3f");
        c->cd(1); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.12); hE->Draw("COLZ TEXT");
        TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.036);
        tl.DrawLatex(0.12,0.94, Form("EARLY — %d GeV (events=%d)", E, usedEvents));
        c->cd(2); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.12); hP->Draw("COLZ TEXT");
        tl.DrawLatex(0.12,0.94, Form("NEAR-PEAK — %d GeV", E));
        c->cd(3); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.12); hL->Draw("COLZ TEXT");
        tl.DrawLatex(0.12,0.94, Form("LATE — %d GeV", E));
        c->SaveAs(Form("lateral_development_realistic_%dGeV.png",E));
        std::printf(" saved lateral_development_realistic_%dGeV.png (events=%d)\n", E, usedEvents);

        // representative near-peak
        TH2D *hRep = (TH2D*)hP->Clone(Form("hREP_%d",E));
        hRep->SetTitle(Form("Representative lateral (near-peak) %d GeV", E));
        TCanvas *ct = new TCanvas(Form("rep_%d",E), "", 600,600);
        hRep->Draw("COLZ TEXT");
        ct->SaveAs(Form("lateral_representative_realistic_%dGeV.png",E));

        // delete heavy objects to free memory (ROOT keeps ownership; deleting helps when re-running)
        delete hE; delete hP; delete hL;
        delete hRep;
        delete c; delete ct;
    } // energies

    Info("plotLateralHeatmaps_safe","Done. energyScale=%.3g perEventNormalize=%s oversample=%d",
         energyScale, perEventNormalize ? "ON":"OFF", oversample);
}

