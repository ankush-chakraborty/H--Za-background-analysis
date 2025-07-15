#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <TLegend.h>     // <-- For TLegend
#include <TStyle.h>      // <-- For gStyle
#include <TROOT.h>       // Optional, for global ROOT settings

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

std::pair<TH1F*, int> histogram(const std::string& filename, const std::string& histname){
    // Open ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
    }

    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        cerr << "Error: Events tree not found!" << endl;
        file->Close();
    }
    
    int c=0;

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
    tree->SetBranchAddress("Particle_energy", E);  // Fixed: was 'energy'
    
    // Histogram for jet pT
    TH1F *hist = new TH1F(histname.c_str(), "; #Delta #phi (Z#gamma, jj) [Radians]; Events", 30, 0, TMath::Pi());
    
    // Loop over events
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        
        vector<PseudoJet> particles;
        vector<pair<int, TLorentzVector>> event_lep;
        vector<TLorentzVector> event_pho;
        
        for (int j = 0; j < nParticles; ++j) {
            // Select only final state particles
            if (status[j] != 1) continue;
            
            //Skip leptons and neutrinos
            int abs_pid = abs(pid[j]);
            
            if (abs_pid == 22) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], E[j]);  // Fixed: was 'energy[j]'
                event_pho.push_back(p4);
            }
            
            if (abs_pid == 11 || abs_pid == 13 || abs_pid == 15 || // e, mu, tau
                abs_pid == 12 || abs_pid == 14 || abs_pid == 16) { // nu_e, nu_mu, nu_tau
                if(abs_pid == 11 || abs_pid == 13) {
                    TLorentzVector p4;
                    p4.SetPxPyPzE(px[j], py[j], pz[j], E[j]);  // Fixed: was 'energy[j]'
                    event_lep.push_back({pid[j], p4});
                    continue;
                } else continue;
            }
            
            // Create the PseudoJet
            fastjet::PseudoJet p(px[j], py[j], pz[j], E[j]);

            // Apply rapidity cut (|η| < 5 is typical for LHC detectors)
            if (std::abs(p.eta()) > 4.7) continue;

            // Optional: skip soft particles
            if (p.perp() < 1.0) continue;

            // Add to particle list
            particles.push_back(p);
        }
        
        // Check if we have enough leptons and photons
        if (event_lep.size() < 2 || event_pho.size() < 1) continue;
        
        // Sort leptons by pT
        sort(event_lep.begin(), event_lep.end(), [](const pair<int, TLorentzVector> &a, const pair<int, TLorentzVector> &b) {
             return a.second.Pt() > b.second.Pt();
        });
        
        // Sort photons by pT
        sort(event_pho.begin(), event_pho.end(), [](const TLorentzVector &a, const TLorentzVector &b) {
             return a.Pt() > b.Pt();
        });
        
        // Check if we have opposite-sign leptons
        if (event_lep[0].first != -event_lep[1].first) continue;
        
        // Form Z+gamma system
        TLorentzVector zg = event_lep[0].second + event_lep[1].second + event_pho[0];

        // Anti-kt jet clustering
        double R = 0.5;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = cs.inclusive_jets(20.0); // minimum jet pT = 20 GeV

        if(jets.size() < 2) continue;
       
        // Sort jets by pT
        sort(jets.begin(), jets.end(), [](const PseudoJet &a, const PseudoJet &b) {
             return a.perp() > b.perp();
        });
        
        PseudoJet j1 = jets[0];
        PseudoJet j2 = jets[1];
        
        if(j2.perp() < 20) continue;
               
        double dy = TMath::Abs(j1.eta() - j2.eta());
        if(dy < 3.5) continue;
              
        double zepp = zg.Eta() - 0.5*(j1.eta() + j2.eta());
        if(zepp > 2.5) continue;
       
        PseudoJet jj = j1 + j2;
        if(jj.m() < 500) continue;
        
        double dphi = TMath::Abs(jj.phi() - zg.Phi());
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        hist->Fill(dphi);
        
        if(dphi < 2.4) continue;
        c++;
    }
    
    return {hist, c};
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
    int count[nhist];
    
    gStyle->SetOptStat(0);  // Disable stats box
    TCanvas *c = new TCanvas("c", "canvas", 800, 600);
    
    for (int i = 0; i < nhist; i++) {
        pair<TH1F*, int> p = histogram(filenames[i], histnames[i]);
        hists[i] = p.first;
        count[i] = p.second;
        if (!hists[i]) {
            cerr << "Error creating histogram " << i << endl;
            continue;
        }
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetLineWidth(2);
        hists[i]->SetTitle("");  // No title
        cout << "Integral: " << hists[i]->Integral() << endl;
        if (hists[i]->Integral() > 0.0) {
            hists[i]->Scale(1.0/hists[i]->Integral());  
        }
    }

    hists[1]->Draw("hist");
    for (int i = 0; i < nhist; i++) {
        if (i == 1) continue;
        hists[i]->Draw("hist same");
    }
    
    // Create and populate legend
    TLegend *legend = new TLegend(0.65, 0.7, 0.88, 0.88);
    for (int i = 0; i < nhist; i++) {
        if (hists[i]) legend->AddEntry(hists[i], legendLabels[i], "l");
    }
    legend->Draw();
    
    c->SaveAs("delta_phi.png");  // Fixed: changed from "delta_eta.png"
    
    std::cout << "ggH : " << count[0] << "\n";
    std::cout << "VBH : " << count[1] << "\n";
    std::cout << "VH : " << count[2] << "\n";
    std::cout << "BKG : " << count[3] << "\n";
    
    // Clean up
    for (int i = 0; i < nhist; i++) {
        if (hists[i]) delete hists[i];
    }
    delete c;
    
    return 0;
}
