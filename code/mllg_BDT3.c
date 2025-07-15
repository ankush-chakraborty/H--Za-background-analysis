#include <TMVA/Reader.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include "TCanvas.h"
#include "TF1.h"      // For defining fit functions
#include "TH1.h"      // For histograms and their Fit() method


void mllg_BDT3() {
    TFile *f = TFile::Open("TMVAOutput_3.root");
    TTree *testTree = (TTree*)f->Get("dataset/TestTree");
    
    float BDTScore;
    int id;
    testTree->SetBranchAddress("BDT", &BDTScore);
    testTree->SetBranchAddress("classID", &id);
    
    std::vector<float> signalScores;
    
    for (Long64_t i = 0; i < testTree->GetEntries(); ++i) {
        testTree->GetEntry(i);
        if (id == 0)
            signalScores.push_back(BDTScore);
    }
    
    std::sort(signalScores.begin(), signalScores.end(), std::greater<float>());
    int cutIndex = signalScores.size() * (0.01 * 70);
    float BDTcut = signalScores[cutIndex];
    std::cout << BDTcut << endl;
    
    TFile *file = TFile::Open("tmva_input1.root");
    TTree *tree = (TTree*)file->Get("tmva_tree");
    
    float deta_jj, dphi_jj, dphi_jj_zg, dR_j_g, pT_j1, pT_j2, pTt, pT_balance, zepp_g, mjj, zepp_h, njets, njets2, pho_pt, l1_pt, l2_pt, pho_eta, z_pt, l_eta;
    int classID, processID;
    float m_llg, theta;
    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  
    reader->AddVariable("deta_jj", &deta_jj);
    reader->AddVariable("dphi_jj", &dphi_jj);
    reader->AddVariable("dphi_jj_zg", &dphi_jj_zg);
    reader->AddVariable("dR_j_g", &dR_j_g);
    reader->AddVariable("pT_j1", &pT_j1);
    reader->AddVariable("pT_j2", &pT_j2);
    reader->AddVariable("pTt", &pTt);
    reader->AddVariable("pT_balance", &pT_balance);
    //reader->AddVariable("zepp_g", &zepp_g);
    reader->AddVariable("zepp_h", &zepp_h);
    //reader->AddVariable("mjj", &mjj);
    //reader->AddVariable("njets", &njets);
    reader->AddVariable("njets2", &njets2);
    //reader->AddVariable("l1_pt", &l1_pt);
    //reader->AddVariable("l2_pt", &l2_pt);
    //reader->AddVariable("pho_pt", &pho_pt);
    //reader->AddVariable("pho_eta", &pho_eta);
    //reader->AddVariable("z_pt", &z_pt);
    //reader->AddVariable("l_eta", &l_eta);
    
    reader->BookMVA("BDT", "dataset/weights/TMVAClassification_BDT.weights.xml");
    
    tree->SetBranchAddress("deta_jj", &deta_jj);
    tree->SetBranchAddress("dphi_jj", &dphi_jj);
    tree->SetBranchAddress("dphi_jj_zg", &dphi_jj_zg);
    tree->SetBranchAddress("dR_j_g", &dR_j_g);
    tree->SetBranchAddress("pT_j1", &pT_j1);
    tree->SetBranchAddress("pT_j2", &pT_j2);
    tree->SetBranchAddress("pTt", &pTt);
    tree->SetBranchAddress("pT_balance", &pT_balance);
    tree->SetBranchAddress("zepp_g", &zepp_g);
    tree->SetBranchAddress("mjj", &mjj);
    tree->SetBranchAddress("zepp_h", &zepp_h);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("njets2", &njets2);
    tree->SetBranchAddress("l1_pt", &l1_pt);
    tree->SetBranchAddress("l2_pt", &l2_pt);
    tree->SetBranchAddress("pho_pt", &pho_pt);
    tree->SetBranchAddress("pho_eta", &pho_eta);
    
    tree->SetBranchAddress("classID", &classID);
    tree->SetBranchAddress("processID", &processID);
    
    tree->SetBranchAddress("m_llg", &m_llg);
    tree->SetBranchAddress("theta", &theta);
    //BDTcut = -999;
    
    int nbins = 25; float xmin = 110, xmax = 150;
    TH1F *h_sig_VBF = new TH1F("h_sig", "Angular Separation of l^{-} and #gamma; #Delta #theta (l^{-}, #gamma) [Radians]; Events", nbins, xmin, xmax);
    TH1F *h_sig_ggH = new TH1F("h_bkg", "Angular Separation of l^{-} and #gamma; #Delta #theta (l^{-}, #gamma) [Radians]; Events", nbins, xmin, xmax);
    TH1F *h_sig_VH = new TH1F("h_bkg", "Angular Separation of l^{-} and #gamma; #Delta #theta (l^{-}, #gamma) [Radians]; Events", nbins, xmin, xmax);
    TH1F *h_bkg = new TH1F("h_bkg", "Angular Separation of l^{-} and #gamma; #Delta #theta (l^{-}, #gamma) [Radians]; Events", nbins, xmin, xmax);
    
    int c_bkg = 0, c_VBF = 0, c_ggH = 0, c_VH = 0;
    int bkg = 0, sig_VBF = 0, sig_ggH = 0, sig_VH = 0;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        
        float bdt_score = reader->EvaluateMVA("BDT");
        if (bdt_score < BDTcut) continue;
        
        bool masscut = (fabs(m_llg - 125) < 1);        
        if (processID == 1) {
        //if (masscut)
            h_sig_VBF->Fill(m_llg);
            c_VBF++;
            if (masscut) sig_VBF++;
        } else if (processID == 2) {
        //if (masscut)
            h_sig_ggH->Fill(m_llg);
            c_ggH++;
            if (masscut) sig_ggH++;
        } else if (processID == 3) {
        //if (masscut)
            h_sig_VH->Fill(m_llg);
            c_VH++;
            if (masscut) sig_VH++;
        } else if (processID == 4) {
        //if (masscut)
            h_bkg->Fill(m_llg);
            c_bkg++;
            if (masscut) bkg++;
        }
    }
    
    std::cout << "Za BKG : " << c_bkg << endl;
    std::cout << "ggH : " << c_ggH << endl;
    std::cout << "VBF : " << c_VBF << endl;
    std::cout << "VH : " << c_VH << endl;
    
    h_sig_VBF->SetLineColor(kRed);
    h_sig_ggH->SetLineColor(kRed);
    h_sig_VH->SetLineColor(kRed);
    h_bkg->SetLineColor(kBlue);
    
    h_sig_VBF->SetFillColor(kRed);
    h_sig_ggH->SetFillColor(kRed);
    h_sig_VH->SetFillColor(kRed);
    h_bkg->SetFillColor(kBlue);
    
    double lum = 300.0;
    double w_ggH = 2.658 * lum / 9995;
    double w_VBF = 0.2604 * lum / 10000;
    double w_VH = 0.1681 * lum / 10000;
    double w_bkg = 5771 * lum / 99602;
    
    h_sig_VBF->Scale(w_VBF);
    h_sig_ggH->Scale(w_ggH);
    h_sig_VH->Scale(w_VH);
    h_bkg->Scale(w_bkg);
    
    /*
    TF1* bkgFit = new TF1("bkgFit", "expo", 110, 150);  // or "pol1", etc.
    h_bkg->Fit(bkgFit, "R");
    
    TH1F* h_sig_total = (TH1F*)h_sig_VBF->Clone("h_sig_total");
    h_sig_total->Add(h_sig_ggH);
    h_sig_total->Add(h_sig_VH);
    
    TF1* sigFit = new TF1("sigFit", "gaus", 110, 150);  // or "pol1", etc.
    h_sig_total->Fit(sigFit, "R");
    
    TF1* totalFit = new TF1("totalFit", "gaus(0) + expo(3)", 110, 150);

    // Set parameters from individual fits
    for (int i = 0; i < 3; ++i)
        totalFit->SetParameter(i, sigFit->GetParameter(i));
    
    for (int i = 0; i < 2; ++i)
        totalFit->SetParameter(i + 3, bkgFit->GetParameter(i));
    
    TH1F* h_sum = (TH1F*) h_sig_total->Clone("h_sum");
    h_sum->Add(h_bkg);
    
    TCanvas* c1 = new TCanvas();
    h_sum->Draw();
    totalFit->SetLineColor(kBlue);
    totalFit->Draw("same");
    
    bkgFit->SetLineColor(kBlue);
    bkgFit->SetLineWidth(2);
    bkgFit->Draw("same");
    sigFit->SetLineColor(kRed);
    sigFit->SetLineWidth(2);
    sigFit->Draw("same");
    */
    
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg->AddEntry(h_sig_VBF, "Signal", "l");
    leg->AddEntry(h_bkg, "Background", "l");
    leg->Draw();
    
    THStack *stack = new THStack("stack", "Cross-section weighted comparison");
    stack->Add(h_bkg);
    //stack->Add(h_sig_ggH);
    //stack->Add(h_sig_VH);
    stack->Add(h_sig_VBF);
    stack->Draw("HIST");
    stack->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");
    stack->GetYaxis()->SetTitle("Events");
    
    gPad->Update();
    c->SaveAs("mllg_3.png");
    
    std::cout << "BKG : " << bkg << endl;
    std::cout << "VBF : " << sig_VBF << endl;
    
    double s = w_ggH * sig_ggH + w_VBF * sig_VBF + w_VH * sig_VH;
    double b = w_bkg * bkg;
    
    std::cout << "Signal : " << s << endl;
    std::cout << "Background : " << b << endl;
    std::cout << "Signal significance : " << s/TMath::Sqrt(s + b) << endl;
}
