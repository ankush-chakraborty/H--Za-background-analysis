#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <TLegend.h>     // <-- For TLegend
#include <TStyle.h>      // <-- For gStyle
#include <TROOT.h>       // Optional, for global ROOT settings

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

TH1F* histogram(const std::string& filename, const std::string& histname){
    // Open ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
    }

    TTree *tree = (TTree*)file->Get("Events");

    // Particle branches
    const int maxP = 10000;
    int nParticles;
    double px[maxP], py[maxP], pz[maxP], E[maxP];
    
    int pid[maxP], status[maxP];
    tree->SetBranchAddress("Particle_pid", pid);
    tree->SetBranchAddress("Particle_status", status);


    tree->SetBranchAddress("Event_numberP", &nParticles);
    tree->SetBranchAddress("Particle_px", px);
    tree->SetBranchAddress("Particle_py", py);
    tree->SetBranchAddress("Particle_pz", pz);
    tree->SetBranchAddress("Particle_energy", E);

    // Histogram for jet pT
    TH1F *hist = new TH1F(histname.c_str(), "; #Delta #eta_{j1, j2}; Events", 30, 0, 10);

    // Loop over events
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        vector<PseudoJet> particles;
        int count = 0;
        for (int j = 0; j < nParticles; ++j) {
            // Select only final state particles
            if (status[j] != 1) continue;

            //Skip leptons and neutrinos
            int abs_pid = abs(pid[j]);
            if (abs_pid == 11 || abs_pid == 13 || abs_pid == 15 || // e, mu, tau
                abs_pid == 12 || abs_pid == 14 || abs_pid == 16)   // nu_e, nu_mu, nu_tau
                {continue;}

            // Create the PseudoJet
            fastjet::PseudoJet p(px[j], py[j], pz[j], E[j]);

            // Apply rapidity cut (|η| < 5 is typical for LHC detectors)
            if (std::abs(p.eta()) > 4.7) continue;

            // Optional: skip soft particles
            if (p.perp() < 1.0) continue;

            // Add to particle list
            particles.push_back(p);
        }

        // Anti-kt jet clustering
        double R = 0.5;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = cs.inclusive_jets(20.0); // minimum jet pT = 20 GeV

        if(jets.size() < 2) continue;
       
        sort(jets.begin(), jets.end(), [](const PseudoJet &a, const PseudoJet &b) {
             return a.perp() > b.perp();
        });
        
        PseudoJet j1 = jets[0];
        PseudoJet j2 = jets[1];
        if(j2.perp() < 20) continue;
               
        double dy = TMath::Abs(j1.eta() - j2.eta());
        hist->Fill(dy);
    }
    return hist;
}

int main() {
    const char* filenames[] = {
        "ggH_showered.root",
        "VBF_showered.root",
        "VH_showered.root",
        "Za_showered.root"
    };
    
    const char* histnames[] = {"hist1", "hist2", "hist3", "hist4"};
    
    const char* legendLabels[] = {
        "ggH",
        "VBF",
        "VH",
        "Z#gamma (Background)"
    };
    
    int colors[] = {kRed, kBlue, kGreen+2, kMagenta};
    
    const int nhist = 4;
    TH1F* hists[nhist];
    
    gStyle->SetOptStat(0);  // Disable stats box
    TCanvas *c = new TCanvas("c", "canvas", 800, 600);
    
    // Load, style, and normalize histograms
    for (int i = 0; i < nhist; i++) {
        hists[i] = histogram(filenames[i], histnames[i]);
        if (!hists[i]) continue;
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetLineWidth(2);
        hists[i]->SetTitle("");  // No title
        hists[i]->Scale(1.0 / hists[i]->Integral());  // Normalize
    }

    hists[2]->Draw("hist");
    for (int i = 0; i < nhist; i++) {
        if (i == 2) continue;
        hists[i]->Draw("hist same");
    }
    
    // Create and populate legend
    TLegend *legend = new TLegend(0.65, 0.7, 0.88, 0.88);
    for (int i = 0; i < nhist; i++) {
        if (hists[i]) legend->AddEntry(hists[i], legendLabels[i], "l");
    }
    legend->Draw();
    
    c->SaveAs("delta_eta.png");
    return 0;
}
