#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <cmath>
#include <fastjet/ClusterSequence.hh>

using namespace std;
using namespace fastjet;

bool is_final_state_hadron(int pid) {
    int apid = abs(pid);
    if (apid == 11 || apid == 12 || apid == 13 || apid == 14 || apid == 15 || apid == 16 || apid == 22)
        return false;
    return true;
}

void create_tmva_input(const string& sig_file, const string& bkg_file1, const string& bkg_file2, const string& bkg_file3, const string& out_file) {
    TFile *fsig = TFile::Open(sig_file.c_str(), "READ");
    if (!fsig || fsig->IsZombie()) {
        cout << "Error: Cannot open signal file " << sig_file << endl;
        return;
    }
    
    TFile *fbkg1 = TFile::Open(bkg_file1.c_str(), "READ");
    if (!fbkg1 || fbkg1->IsZombie()) {
        cout << "Error: Cannot open background file " << bkg_file1 << endl;
        fsig->Close();
        return;
    }
    
    TFile *fbkg2 = TFile::Open(bkg_file2.c_str(), "READ");
    if (!fbkg2 || fbkg2->IsZombie()) {
        cout << "Error: Cannot open background file " << bkg_file2 << endl;
        fsig->Close();
        return;
    }
    
    TFile *fbkg3 = TFile::Open(bkg_file3.c_str(), "READ");
    if (!fbkg3 || fbkg3->IsZombie()) {
        cout << "Error: Cannot open background file " << bkg_file3 << endl;
        fsig->Close();
        return;
    }

    TTree *tsig = (TTree*)fsig->Get("Events");
    TTree *tbkg1 = (TTree*)fbkg1->Get("Events");
    TTree *tbkg2 = (TTree*)fbkg2->Get("Events");
    TTree *tbkg3 = (TTree*)fbkg3->Get("Events");
    
    if (!tsig || !tbkg1 || !tbkg2 || !tbkg3) {
        fsig->Close();
        fbkg1->Close();
        fbkg2->Close();
        fbkg3->Close();
        return;
    }
    
    const int MAX_PARTICLES = 5000;
    Int_t nparticles;
    Double_t px[MAX_PARTICLES], py[MAX_PARTICLES], pz[MAX_PARTICLES], E[MAX_PARTICLES];
    Int_t pid[MAX_PARTICLES], status[MAX_PARTICLES];
    
    auto set_branches = [&](TTree* tree) {
        if (tree->GetBranch("Event_numberP")) tree->SetBranchAddress("Event_numberP", &nparticles);
        else { cout << "Warning: Branch 'Event_numberP' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_px")) tree->SetBranchAddress("Particle_px", px);
        else { cout << "Warning: Branch 'Particle_px' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_py")) tree->SetBranchAddress("Particle_py", py);
        else { cout << "Warning: Branch 'Particle_py' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_pz")) tree->SetBranchAddress("Particle_pz", pz);
        else { cout << "Warning: Branch 'Particle_pz' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_energy")) tree->SetBranchAddress("Particle_energy", E);
        else { cout << "Warning: Branch 'Particle_energy' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_pid")) tree->SetBranchAddress("Particle_pid", pid);
        else { cout << "Warning: Branch 'Particle_pid' not found" << endl; return false; }
        
        if (tree->GetBranch("Particle_status")) tree->SetBranchAddress("Particle_status", status);
        else { cout << "Warning: Branch 'Particle_status' not found" << endl; return false; }
        
        return true;
    };
    
    if (!set_branches(tsig) || !set_branches(tbkg1) || !set_branches(tbkg2) || !set_branches(tbkg3)) {
        fsig->Close();
        fbkg1->Close();
        fbkg2->Close();
        fbkg3->Close();
        return;
    }
    
    TFile *fout = new TFile(out_file.c_str(), "RECREATE");
    TTree *tout = new TTree("tmva_tree", "TMVA input tree");
    
    float deta_jj, dphi_jj, dphi_jj_zg, dR_j_g, pT_j1, pT_j2, pTt, pT_balance, zepp_g, zepp_h, mjj, weight, pho_pt, l1_pt, l2_pt, pho_eta, z_pt, l_eta, var1, var2, var3;
    int njets, classID, processID, njets2;
    float m_llg, theta;
    
    tout->Branch("deta_jj", &deta_jj);
    tout->Branch("dphi_jj", &dphi_jj);
    tout->Branch("dphi_jj_zg", &dphi_jj_zg);
    tout->Branch("dR_j_g", &dR_j_g);
    tout->Branch("pT_j1", &pT_j1);
    tout->Branch("pT_j2", &pT_j2);
    tout->Branch("pTt", &pTt);
    tout->Branch("pT_balance", &pT_balance);
    tout->Branch("zepp_g", &zepp_g);
    tout->Branch("zepp_h", &zepp_h);
    tout->Branch("classID", &classID);
    tout->Branch("processID", &processID);
    
    tout->Branch("mjj", &mjj);
    tout->Branch("njets", &njets);
    tout->Branch("njets2", &njets2);
    
    tout->Branch("l1_pt", &l1_pt);
    tout->Branch("l2_pt", &l2_pt);
    tout->Branch("pho_pt", &pho_pt);
    tout->Branch("pho_eta", &pho_eta);
    tout->Branch("z_pt", &z_pt);
    tout->Branch("l_eta", &l_eta);
    tout->Branch("var1", &var1);
    tout->Branch("var2", &var2);
    tout->Branch("var3", &var3);
    
    tout->Branch("m_llg", &m_llg);
    tout->Branch("theta", &theta);
    
    tout->Branch("weight", &weight);
    
    auto process_tree = [&](TTree *tree, int label, int procID) {
        Long64_t nentries = tree->GetEntries();
        cout << "Processing " << nentries << " events for process " << procID << endl;
        
        JetDefinition jet_def(antikt_algorithm, 0.4);

        for (Long64_t i = 0; i < nentries; ++i) {
            if (i % 1000 == 0) cout << "Processing event " << i << "/" << nentries << endl;
            
            tree->GetEntry(i);
            
            // Bounds checking for nparticles
            if (nparticles <= 0 || nparticles > MAX_PARTICLES) {
                cout << "Warning: Invalid nparticles = " << nparticles << " in event " << i << endl;
                continue;
            }
            
            vector<PseudoJet> fj_particles;
            vector<pair<int, TLorentzVector>> leptons;
            vector<TLorentzVector> photons;
            fj_particles.reserve(nparticles); // Reserve space to avoid reallocations

            for (int j = 0; j < nparticles; ++j) {
                if (status[j] != 1) continue;
                
                if (!is_final_state_hadron(pid[j])) {
                    if (abs(pid[j]) == 22) {
                        TLorentzVector p;
                        p.SetPxPyPzE(px[j], py[j], pz[j], E[j]);
                        photons.push_back(p);
                        continue;
                    }
                    
                    if (abs(pid[j]) == 11 || abs(pid[j]) == 13) {
                        TLorentzVector p;
                        p.SetPxPyPzE(px[j], py[j], pz[j], E[j]);
                        leptons.push_back({pid[j], p});
                        continue;
                    }
                }
                
                if (!isfinite(px[j]) || !isfinite(py[j]) || !isfinite(pz[j]) || !isfinite(E[j])) {
                    cout << "Warning: Invalid momentum values in event " << i << ", particle " << j << endl;
                    continue;
                }
                
                PseudoJet p(px[j], py[j], pz[j], E[j]);              
                fj_particles.push_back(p);
            }
            
            if (fj_particles.size() < 2) continue;
            
            ClusterSequence cs(fj_particles, jet_def);
            vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(30.0));
            
            sort(leptons.begin(), leptons.end(), [](const pair<int, TLorentzVector> &a, const pair<int, TLorentzVector> &b) {
                return a.second.Pt() > b.second.Pt();
            });
            
            sort(photons.begin(), photons.end(), [](const TLorentzVector &a, const TLorentzVector &b) {
                return a.Pt() > b.Pt();
            });
            
            if (jets.size() < 2) continue;
            
            TLorentzVector j1, j2;
            j1.SetPxPyPzE(jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].e());
            j2.SetPxPyPzE(jets[1].px(), jets[1].py(), jets[1].pz(), jets[1].e());
            
            for(int i = 1; i < leptons.size(); i++) {
                if (leptons[0].first != -leptons[i].first) continue;
                leptons[1] = leptons[i];
                break;
            }
            
            TLorentzVector l1, l2;
            if (leptons[0].first == 11) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else if (leptons[1].first == 11) {l1 = leptons[1].second; l2 = leptons[0].second;}
            if (leptons[0].first == 13) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else if (leptons[1].first == 13) {l1 = leptons[1].second; l2 = leptons[0].second;}
            if (leptons[0].first == 15) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else if (leptons[1].first == 15) {l1 = leptons[1].second; l2 = leptons[0].second;}
            TLorentzVector gamma = photons[0];
            
            if (fabs(j1.Eta()) > 4.7 || fabs(j2.Eta()) > 4.7) continue;
            
            if (j1.DeltaR(l1) < 0.4 || j1.DeltaR(l2) < 0.4 || j1.DeltaR(gamma) < 0.4) continue;
            if (j2.DeltaR(l1) < 0.4 || j2.DeltaR(l2) < 0.4 || j2.DeltaR(gamma) < 0.4) continue;
            
            deta_jj = fabs(j1.Eta() - j2.Eta());
            dphi_jj = fabs(j1.Phi() - j2.Phi());
            if (dphi_jj > TMath::Pi()) dphi_jj = 2*TMath::Pi() - dphi_jj;
            dphi_jj_zg = fabs((j1 + j2).Phi() - (l1 + l2 + gamma).Phi());
            if (dphi_jj_zg > TMath::Pi()) dphi_jj_zg = 2*TMath::Pi() - dphi_jj_zg;
            dR_j_g = min(j1.DeltaR(gamma), j2.DeltaR(gamma));
            pT_j1 = j1.Pt();
            pT_j2 = j2.Pt();
            
            double pT_H = (l1 + l2 + gamma).Pt();
            double pT_zg_cross = (l1 + l2).Px() * gamma.Py() - (l1 + l2).Py() * gamma.Px();
            pTt = 2/pT_H * fabs(pT_zg_cross);
            
            double a = (j1 + gamma).Pt();
            double b = j1.Px() * gamma.Py() - j1.Py() * gamma.Px();
            var1 = 2/a * abs(b);
            
            a = (j2 + gamma).Pt();
            b = j2.Px() * gamma.Py() - j2.Py() * gamma.Px();
            var2 = 2/a * abs(b);
            
            a = (j1 + j2 + gamma).Pt();
            b = (j1 + j2).Px() * gamma.Py() - (j1 + j2).Py() * gamma.Px();
            var3 = 2/a * abs(b);
            
            TLorentzVector p_sum = l1 + l2 + gamma + j1 + j2;
            double pT_vec_sum = fabs(p_sum.Pt());
            double pT_sum = fabs((l1 + l2).Pt() + gamma.Pt() + j1.Pt() + j2.Pt());
            pT_balance = pT_vec_sum/pT_sum;
            
            zepp_g = fabs(gamma.Eta() - 0.5*(j1.Eta() + j2.Eta()));
            zepp_h = fabs((l1 + l2 + gamma).Eta() - 0.5*(j1.Eta() + j2.Eta()));
            
            njets = jets.size();
            njets2 = 0;
            if (jets.size() == 2) njets2 = 0;
            else {
                for (int i = 2; i < jets.size(); i++) {
                    if (jets[i].eta() < min(j1.Eta(), j2.Eta()) || jets[i].eta() > max(j1.Eta(), j2.Eta())) continue;
                    njets2++;
                }
            }
                        
            mjj = (j1 + j2).M();
            m_llg = (l1 + l2 + gamma).M();
            pho_pt = gamma.Pt();
            l1_pt = l1.Pt();
            l2_pt = l2.Pt();
            pho_eta = fabs(gamma.Eta());
            l_eta = max(fabs(l1.Eta()), fabs(l2.Eta()));
            z_pt = (l1 + l2).Pt();
            
            TLorentzVector z = l1 + l2;
            TLorentzVector h = z + gamma;
            TVector3 beta1 = z.BoostVector();
            TVector3 beta2 = h.BoostVector();
            
            l1.Boost(-beta1);
            l2.Boost(-beta1);
            z.Boost(-beta2);
            
            theta = l1.Angle(z.Vect());
            
            classID = label;
            processID = procID;
            
            if (processID == 1) weight = 2.658 * 300 / 9995;
            if (processID == 2) weight = 0.2604 * 300 / 10000;
            if (processID == 3) weight = 0.1681 * 300 / 10000;
            if (processID == 4) weight = 5771 * 300 / 99602;
            tout->Fill();
        }
    };

    process_tree(tsig, 1, 1);
    //process_tree(tbkg1, 1, 2);
    //process_tree(tbkg2, 1, 3);
    process_tree(tbkg3, 0, 4);

    fout->cd();
    tout->Write();
    fout->Close();
    
    fsig->Close();
    fbkg1->Close();
    fbkg2->Close();
    fbkg3->Close();
    
    cout << "TMVA input file saved to: " << out_file << endl;
}

int main() {
    create_tmva_input("VBF_showered.root",
                      "ggH_showered.root",
                      "VH_showered.root",
                      "Za_showered.root",
                      "tmva_input1.root");
    return 0;
}
