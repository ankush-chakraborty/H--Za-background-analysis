#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

struct Event {
    vector<PseudoJet> jets;
    vector<pair<int, TLorentzVector>> leptons;
    vector<TLorentzVector> photon;
    int number;
};

double dR(const PseudoJet& j, const TLorentzVector& p) {
    double deta = j.eta() - p.Eta();
    double dphi = j.phi() - p.Phi();
    double r = TMath::Sqrt(deta*deta + dphi*dphi);
    return r;
}

double et(const PseudoJet& jet) {
    return TMath::Sqrt(jet.perp()*jet.perp() + jet.m()*jet.m());
}

double Et(const TLorentzVector& particle) {
    return TMath::Sqrt(particle.Pt()*particle.Pt() + particle.M()*particle.M());
}

bool is_lepton_tag(const Event& event) {
    if (event.leptons.size() < 3) return false;
    pair<int, TLorentzVector> lepton3 = event.leptons[2];
    
    if (lepton3.second.Pt() < 7 && abs(lepton3.first) == 11) return false;
    if (lepton3.second.Pt() < 5 && abs(lepton3.first) == 13) return false;
    return true;
}

bool is_dijet_tag(const Event& event) {
    vector<PseudoJet> jets = event.jets;
    if (jets.size() < 2) return false;
    
    sort(jets.begin(), jets.end(), [](const PseudoJet &a, const PseudoJet &b) {
        return et(a) > et(b);
    });
    
    PseudoJet j1 = jets[0];
    PseudoJet j2 = jets[1];
    PseudoJet jj = j1 + j2;
    
    TLorentzVector l1 = event.leptons[0].second;
    TLorentzVector l2 = event.leptons[0].second;
    TLorentzVector p1 = event.photon[0];
    TLorentzVector zg = l1 + l2 + p1;
    
    if (dR(j1, l1) < 0.4 || dR(j1, l2) < 0.4 || dR(j1, p1) < 0.4 ||
        dR(j2, l1) < 0.4 || dR(j2, l2) < 0.4 || dR(j2, p1) < 0.4) return false;
        
    if (et(j2) < 30) return false;
    
    if (j1.rapidity() > 4.7 || j2.rapidity() > 4.7) return false;
    
    if (abs(j1.rapidity() - j2.rapidity()) < 3.5) return false;
    
    if (zg.Rapidity() - 0.5*(j1.rapidity() + j2.rapidity()) > 2.5) return false;
    
    if (jj.m() < 500) return false;
    
    if (TMath::Abs(zg.Phi() - jj.phi()) < 2.4) return false;
    
    return true;    
}

bool is_boosted_tag(const Event& event){
    TLorentzVector l1 = event.leptons[0].second;
    TLorentzVector l2 = event.leptons[1].second;
    TLorentzVector p1 = event.photon[0];
    TLorentzVector zg = l1 + l2 + p1;
    
    if (zg.Pt() < 60) return false;
    
    return true;
}

tuple<vector<int>, vector<int>, vector<int>, vector<int>> categorize(const string& filename) {

    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
    }

    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        cerr << "Error: Events tree not found!" << endl;
        file->Close();
    }
    
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
    
    vector<int> lepton_tag;
    vector<int> dijet_tag;
    vector<int> boosted_tag;
    vector<int> untagged;
    
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        Event event;
        tree->GetEntry(i);
        
        vector<PseudoJet> particles;
        vector<pair<int, TLorentzVector>> event_lep;
        vector<TLorentzVector> event_pho;
        
        for (int j = 0; j < nParticles; ++j) {
            if (status[j] != 1) continue;
            
            int abs_pid = abs(pid[j]);
            
            if (abs_pid == 22) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], E[j]);
                event_pho.push_back(p4);
            }
            
            if (abs_pid == 11 || abs_pid == 13 || abs_pid == 15 || 
                abs_pid == 12 || abs_pid == 14 || abs_pid == 16) {
                if(abs_pid == 11 || abs_pid == 13) {
                    TLorentzVector p4;
                    p4.SetPxPyPzE(px[j], py[j], pz[j], E[j]);
                    event_lep.push_back({pid[j], p4});
                    continue;
                } else continue;
            }
            
            PseudoJet p(px[j], py[j], pz[j], E[j]);
            if (abs(p.eta()) > 4.7) continue;
            
            if (p.perp() < 1.0) continue;
            
            particles.push_back(p);
        }
        
        sort(event_lep.begin(), event_lep.end(), [](const pair<int, TLorentzVector> &a, const pair<int, TLorentzVector> &b) {
             return a.second.Pt() > b.second.Pt();
        });
        
        sort(event_pho.begin(), event_pho.end(), [](const TLorentzVector &a, const TLorentzVector &b) {
             return a.Pt() > b.Pt();
        });
        
        bool found_pair = false;
        for (int i = 1; i < event_lep.size(); i++) {
            if (event_lep[0].first == -event_lep[i].first) {
                swap(event_lep[1], event_lep[i]);
                found_pair = true;
                break;
            }
        }
        if (!found_pair) continue;
        
        if(event_lep.size() > 0) event.leptons.push_back(event_lep[0]);
        if(event_lep.size() > 1) event.leptons.push_back(event_lep[1]);
        if(event_lep.size() > 2) event.leptons.push_back(event_lep[2]);
        if(event_pho.size() > 0) event.photon.push_back(event_pho[0]);
        
        double R = 0.4;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = cs.inclusive_jets(20.0);
       
        sort(jets.begin(), jets.end(), [](const PseudoJet &a, const PseudoJet &b) {
             return et(a) > et(b);
        });
        
        if(jets.size() > 0) event.jets.push_back(jets[0]);
        if(jets.size() > 1) event.jets.push_back(jets[1]);
        event.number = i;
        
        if(is_lepton_tag(event)) {
            lepton_tag.push_back(event.number);
        } else if (is_dijet_tag(event)) {
            dijet_tag.push_back(event.number);
        } else if (is_boosted_tag(event)) {
            boosted_tag.push_back(event.number);
        } else {
            untagged.push_back(event.number);
        }
        
    }
    
    return make_tuple(lepton_tag, dijet_tag, boosted_tag, untagged);
}

int main() {
    const char* filenames[] = {"ggH_showered.root", "VBF_showered.root", "VH_showered.root", "Za_showered.root"};
    
    vector<int> lepton_tag, dijet_tag, boosted_tag, untagged;
    for (int i = 0; i < 4; i++) {    
        tuple<vector<int>, vector<int>, vector<int>, vector<int>> list = categorize(filenames[i]);
        lepton_tag = get<0>(list);
        dijet_tag = get<1>(list);
        boosted_tag = get<2>(list);
        untagged = get<3>(list);
        
        cout << filenames[i] << endl;
        cout << "Lepton-tag : " << lepton_tag.size() << endl;
        cout << "Dijet-tag : " << dijet_tag.size() << endl;
        cout << "Boosted-tag : " << boosted_tag.size() << endl;
        cout << "Untagged : " << untagged.size() << endl;
    }
    
    return 0;
}
