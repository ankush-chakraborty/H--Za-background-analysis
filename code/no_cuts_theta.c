void no_cuts_1() {
    int nbins = 25; float xmin = 0, xmax = TMath::Pi();
    int bkg = 0, sig_VBF = 0, sig_ggH = 0, sig_VH = 0;
    TFile *file1 = TFile::Open("ggH_showered.root");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree *tree1 = (TTree*)file1->Get("Events");
    if (!tree1) {
        std::cerr << "Error: TTree 'Events' not found." << std::endl;
        return;
    }
    
    TFile *file2 = TFile::Open("VBF_showered.root");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree *tree2 = (TTree*)file2->Get("Events");
    if (!tree2) {
        std::cerr << "Error: TTree 'Events' not found." << std::endl;
        return;
    }
    
    TFile *file3 = TFile::Open("VH_showered.root");
    if (!file3 || file3->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree *tree3 = (TTree*)file3->Get("Events");
    if (!tree3) {
        std::cerr << "Error: TTree 'Events' not found." << std::endl;
        return;
    }
    
    TFile *file4 = TFile::Open("Za_showered.root");
    if (!file4 || file4->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree *tree4 = (TTree*)file4->Get("Events");
    if (!tree4) {
        std::cerr << "Error: TTree 'Events' not found." << std::endl;
        return;
    }

    const int MAXP = 5000;
    Int_t nParticles;
    Int_t pid[MAXP], status[MAXP];
    Double_t px[MAXP], py[MAXP], pz[MAXP], energy[MAXP];

    tree1->SetBranchAddress("Event_numberP", &nParticles); // double-check this name
    tree1->SetBranchAddress("Particle_pid", pid);
    tree1->SetBranchAddress("Particle_status", status);
    tree1->SetBranchAddress("Particle_px", px);
    tree1->SetBranchAddress("Particle_py", py);
    tree1->SetBranchAddress("Particle_pz", pz);
    tree1->SetBranchAddress("Particle_energy", energy);

    TH1F *h1 = new TH1F("h", "Angular Separation; #theta [Radians]; Events", nbins, xmin, xmax);

    Long64_t nentries = tree1->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree1->GetEntry(i);

        std::vector<std::pair<float, TLorentzVector>> leptons;
        std::vector<std::pair<float, TLorentzVector>> photons;

        for (int j = 0; j < nParticles && j < MAXP; ++j) {
            if (status[j] == 1 && (abs(pid[j]) == 11 || abs(pid[j]) == 13)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                leptons.push_back({pid[j], p4});
            }
            
            if (status[j] == 1 && (abs(pid[j]) == 22)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                photons.push_back({pid[j], p4});
            }
        }

        std::sort(leptons.begin(), leptons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        std::sort(photons.begin(), photons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        
        bool found_pair = false;
        for (int i = 1; i < leptons.size(); i++) {
            if (leptons[0].first == -leptons[i].first) {
                swap(leptons[1], leptons[i]);
                found_pair = true;
                break;
            }
        }
        if (!found_pair) continue;

        if (leptons.size() >= 2 && photons.size() > 1) {
            TLorentzVector l1, l2;
            if(leptons[0].first > 0) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else {l2 = leptons[0].second; l1 = leptons[1].second;}
            TLorentzVector gamma = photons[0].second;
            
            TLorentzVector z = l1 + l2;
            TLorentzVector h = z + gamma;
            TVector3 b1 = z.BoostVector();
            TVector3 b2 = h.BoostVector();
            l1.Boost(-b1);
            z.Boost(-b2);
            
            double theta = l1.Angle(z.Vect());
            if(fabs(h.M() - 125) < 2) {sig_ggH++;
            h1->Fill(theta);}
        }
    }
    
    tree2->SetBranchAddress("Event_numberP", &nParticles); // double-check this name
    tree2->SetBranchAddress("Particle_pid", pid);
    tree2->SetBranchAddress("Particle_status", status);
    tree2->SetBranchAddress("Particle_px", px);
    tree2->SetBranchAddress("Particle_py", py);
    tree2->SetBranchAddress("Particle_pz", pz);
    tree2->SetBranchAddress("Particle_energy", energy);

    TH1F *h2 = new TH1F("h", "Angular Separation; #theta [Radians]; Events", nbins, xmin, xmax);

    nentries = tree2->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree2->GetEntry(i);

        std::vector<std::pair<float, TLorentzVector>> leptons;
        std::vector<std::pair<float, TLorentzVector>> photons;

        for (int j = 0; j < nParticles && j < MAXP; ++j) {
            if (status[j] == 1 && (abs(pid[j]) == 11 || abs(pid[j]) == 13)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                leptons.push_back({pid[j], p4});
            }
            
            if (status[j] == 1 && (abs(pid[j]) == 22)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                photons.push_back({pid[j], p4});
            }
        }

        std::sort(leptons.begin(), leptons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        std::sort(photons.begin(), photons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        
        bool found_pair = false;
        for (int i = 1; i < leptons.size(); i++) {
            if (leptons[0].first == -leptons[i].first) {
                swap(leptons[1], leptons[i]);
                found_pair = true;
                break;
            }
        }
        if (!found_pair) continue;

        if (leptons.size() >= 2 && photons.size() > 1) {
            TLorentzVector l1, l2;
            if(leptons[0].first > 0) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else {l2 = leptons[0].second; l1 = leptons[1].second;}
            TLorentzVector gamma = photons[0].second;
            
            TLorentzVector z = l1 + l2;
            TLorentzVector h = z + gamma;
            TVector3 b1 = z.BoostVector();
            TVector3 b2 = h.BoostVector();
            l1.Boost(-b1);
            z.Boost(-b2);
            
            double theta = l1.Angle(z.Vect());
            if(fabs(h.M() - 125) < 2) {sig_VBF++;
            h2->Fill(theta);}
        }
    }
    
    tree3->SetBranchAddress("Event_numberP", &nParticles); // double-check this name
    tree3->SetBranchAddress("Particle_pid", pid);
    tree3->SetBranchAddress("Particle_status", status);
    tree3->SetBranchAddress("Particle_px", px);
    tree3->SetBranchAddress("Particle_py", py);
    tree3->SetBranchAddress("Particle_pz", pz);
    tree3->SetBranchAddress("Particle_energy", energy);

    TH1F *h3 = new TH1F("h", "Angular Separation; #theta [Radians]; Events", nbins, xmin, xmax);

    nentries = tree3->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree3->GetEntry(i);

        std::vector<std::pair<float, TLorentzVector>> leptons;
        std::vector<std::pair<float, TLorentzVector>> photons;

        for (int j = 0; j < nParticles && j < MAXP; ++j) {
            if (status[j] == 1 && (abs(pid[j]) == 11 || abs(pid[j]) == 13)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                leptons.push_back({pid[j], p4});
            }
            
            if (status[j] == 1 && (abs(pid[j]) == 22)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                photons.push_back({pid[j], p4});
            }
        }

        std::sort(leptons.begin(), leptons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        std::sort(photons.begin(), photons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        
        bool found_pair = false;
        for (int i = 1; i < leptons.size(); i++) {
            if (leptons[0].first == -leptons[i].first) {
                swap(leptons[1], leptons[i]);
                found_pair = true;
                break;
            }
        }
        if (!found_pair) continue;

        if (leptons.size() >= 2 && photons.size() > 1) {
            TLorentzVector l1, l2;
            if(leptons[0].first > 0) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else {l2 = leptons[0].second; l1 = leptons[1].second;}
            TLorentzVector gamma = photons[0].second;
            
            TLorentzVector z = l1 + l2;
            TLorentzVector h = z + gamma;
            TVector3 b1 = z.BoostVector();
            TVector3 b2 = h.BoostVector();
            l1.Boost(-b1);
            z.Boost(-b2);
            
            double theta = l1.Angle(z.Vect());
            if(fabs(h.M() - 125) < 2) {sig_VH++;
            h3->Fill(theta);}
        }
    }
    
    tree4->SetBranchAddress("Event_numberP", &nParticles); // double-check this name
    tree4->SetBranchAddress("Particle_pid", pid);
    tree4->SetBranchAddress("Particle_status", status);
    tree4->SetBranchAddress("Particle_px", px);
    tree4->SetBranchAddress("Particle_py", py);
    tree4->SetBranchAddress("Particle_pz", pz);
    tree4->SetBranchAddress("Particle_energy", energy);

    TH1F *h4 = new TH1F("h", "Angular Separation; #theta [Radians]; Events", nbins, xmin, xmax);

    nentries = tree4->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree4->GetEntry(i);

        std::vector<std::pair<float, TLorentzVector>> leptons;
        std::vector<std::pair<float, TLorentzVector>> photons;

        for (int j = 0; j < nParticles && j < MAXP; ++j) {
            if (status[j] == 1 && (abs(pid[j]) == 11 || abs(pid[j]) == 13)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                leptons.push_back({pid[j], p4});
            }
            
            if (status[j] == 1 && (abs(pid[j]) == 22)) {
                TLorentzVector p4;
                p4.SetPxPyPzE(px[j], py[j], pz[j], energy[j]);
                photons.push_back({pid[j], p4});
            }
        }

        std::sort(leptons.begin(), leptons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        std::sort(photons.begin(), photons.end(), [](auto &a, auto &b) { return a.second.Pt() > b.second.Pt(); });
        
        //for (int i = 1; i < leptons.size(); i++){
          //  if (leptons[0].first != -leptons[i].first) continue;
            //leptons[1] = leptons[i];
        //}

        if (leptons.size() >= 2 && photons.size() >= 1) {
            TLorentzVector l1, l2;
            if(leptons[0].first > 0) {l1 = leptons[0].second; l2 = leptons[1].second;}
            else {l2 = leptons[0].second; l1 = leptons[1].second;}
            TLorentzVector gamma = photons[0].second;
            
            TLorentzVector z = l1 + l2;
            TLorentzVector h = z + gamma;
            TVector3 b1 = z.BoostVector();
            TVector3 b2 = h.BoostVector();
            l1.Boost(-b1);
            z.Boost(-b2);
            
            double theta = l1.Angle(z.Vect());
            if(fabs(h.M() - 125) < 2) {bkg++;
            h4->Fill(theta);}
        }
    }
    
    double lum = 300.0;
    double w_ggH = 2.658 * lum / 9995;
    double w_VBF = 0.2604 * lum / 10000;
    double w_VH = 0.1681 * lum / 10000;
    double w_za = 5771 * lum / 99602;
    
    //w_ggH = 1;
    //w_za = 3;
    
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    THStack *stack = new THStack("stack", "Cross-section weighted comparison");
    
    h1->Scale(w_ggH);
    h2->Scale(w_VBF);
    h3->Scale(w_VH);
    h4->Scale(w_za);
    
    h1->SetLineColor(kRed);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kRed);
    h4->SetLineColor(kBlue);
    
    h1->SetFillColor(kRed);
    h2->SetFillColor(kRed);
    h3->SetFillColor(kRed);
    h4->SetFillColor(kBlue);
    
    stack->Add(h4);
    //stack->Add(h3);
    stack->Add(h2);
    //stack->Add(h1);
    stack->Draw("HIST");
    
    double s = w_ggH * sig_ggH + w_VBF * sig_VBF + w_VH * sig_VH;
    double b = w_za * bkg;
    
    std::cout << "Signal : " << s << endl;
    std::cout << "Background : " << b << endl;
    std::cout << "Signal significance : " << s/TMath::Sqrt(s + b) << endl;
    c->SaveAs("p1223.png");
}
