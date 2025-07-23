#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TBranch.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

bool is_final_state_hadron(int pid) {
    int apid = abs(pid);
    if (apid == 11 || apid == 12 || apid == 13 || apid == 14 || apid == 15 || apid == 16 || apid == 22)
        return false;
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_file>\n";
        return 1;
    }

    std::string inputFilename = argv[1];
    TFile *inputFile = new TFile(inputFilename.c_str());
    TTree *inputTree = (TTree*) inputFile->Get("Events");  // Replace with actual TTree name

    // Input variables
    const int MAX_PARTICLES = 5000;
    Int_t nparticles;
    Double_t px[MAX_PARTICLES], py[MAX_PARTICLES], pz[MAX_PARTICLES], E[MAX_PARTICLES];
    Int_t pid[MAX_PARTICLES], status[MAX_PARTICLES];

    inputTree->SetBranchAddress("Particle_px", &px);
    inputTree->SetBranchAddress("Particle_py", &py);
    inputTree->SetBranchAddress("Particle_pz", &pz);
    inputTree->SetBranchAddress("Particle_energy",  &E);
    inputTree->SetBranchAddress("Particle_pid", &pid);
    inputTree->SetBranchAddress("Event_numberP", &nparticles);
    inputTree->SetBranchAddress("Particle_status", &status);

    // Output file and tree
    TFile *outputFile = new TFile("output.root", "RECREATE");
    TTree *outTree = new TTree("events", "Reconstructed llgamma and jets");

    // Output variables
    vector<int> out_pid;
    vector<float> out_px, out_py, out_pz, out_E;

    // Create branches
    outTree->Branch("pid", &out_pid);
    outTree->Branch("px", &out_px);
    outTree->Branch("py", &out_py);
    outTree->Branch("pz", &out_pz);
    outTree->Branch("energy", &out_E);

    Long64_t nentries = inputTree->GetEntries();
    JetDefinition jet_def(antikt_algorithm, 0.4);

    for (Long64_t i = 0; i < nentries; ++i) {
        inputTree->GetEntry(i);

        if (nparticles <= 0 || nparticles > MAX_PARTICLES) {
            cout << "Warning: Invalid nparticles = " << nparticles << " in event " << i << endl;
            continue;
        }

        vector<PseudoJet> fj_particles;
        vector<pair<int, TLorentzVector>> leptons;
        vector<TLorentzVector> photons;

        fj_particles.reserve(nparticles);

        for (int j = 0; j < nparticles; ++j) {
            if (status[j] != 1) continue;

            if (!isfinite(px[j]) || !isfinite(py[j]) || !isfinite(pz[j]) || !isfinite(E[j])) {
                cout << "Warning: Invalid momentum values in event " << i << ", particle " << j << endl;
                continue;
            }

            if (!is_final_state_hadron(pid[j])) {
                TLorentzVector p;
                p.SetPxPyPzE(px[j], py[j], pz[j], E[j]);

                if (abs(pid[j]) == 22) {
                    photons.push_back(p);
                } else if (abs(pid[j]) == 11 || abs(pid[j]) == 13) {
                    leptons.push_back({pid[j], p});
                }
            } else {
                fj_particles.emplace_back(px[j], py[j], pz[j], E[j]);
            }
        }

        vector<PseudoJet> jets;
        if (fj_particles.size() >= 2) {
            ClusterSequence cs(fj_particles, jet_def);
            jets = sorted_by_pt(cs.inclusive_jets(30.0));
        }

        // Clear output vectors
        out_pid.clear();
        out_px.clear();
        out_py.clear();
        out_pz.clear();
        out_E.clear();

        // Fill leptons
        for (auto &lep : leptons) {
            out_pid.push_back(lep.first);
            out_px.push_back(lep.second.Px());
            out_py.push_back(lep.second.Py());
            out_pz.push_back(lep.second.Pz());
            out_E.push_back(lep.second.E());
        }

        // Fill photons
        for (auto &pho : photons) {
            out_pid.push_back(22);
            out_px.push_back(pho.Px());
            out_py.push_back(pho.Py());
            out_pz.push_back(pho.Pz());
            out_E.push_back(pho.E());
        }
        
        for (auto &jet : jets) {
            out_pid.push_back(999);
            out_px.push_back(jet.px());
            out_py.push_back(jet.py());
            out_pz.push_back(jet.pz());
            out_E.push_back(jet.e());
        }

        // Fill tree once per event
        outTree->Fill();
    }

    // Save output
    outputFile->cd();
    outTree->Write();
    outputFile->Close();
    inputFile->Close();

    cout << "Done. Output written to reco_output.root" << endl;

    return 0;
}

