void plot() {
    TMVA::Tools::Instance();

    // Open files
    TFile *file1 = TFile::Open("TMVAOutput_1.root");
    TFile *file2 = TFile::Open("TMVAOutput_2.root");
    TFile *file3 = TFile::Open("TMVAOutput_3.root");

    // Get ROC curves
    TGraph* roc1 = (TGraph*)file1->Get("dataset/Method_BDT/BDT/MVA_BDT_rejBvsS");
    TGraph* roc2 = (TGraph*)file2->Get("dataset/Method_BDT/BDT/MVA_BDT_rejBvsS");
    TGraph* roc3 = (TGraph*)file3->Get("dataset/Method_BDT/BDT/MVA_BDT_rejBvsS");

    // Set graph styles
    roc1->SetLineColor(kRed);
    roc1->SetLineWidth(2);
    roc2->SetLineColor(kMagenta);
    roc2->SetLineWidth(2);
    roc3->SetLineColor(kBlue);
    roc3->SetLineWidth(2);

    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "ROC Curves", 800, 600);

    // Draw dummy histogram to provide axes
    TH1F *frame = new TH1F("frame", "ROC Comparison;Signal Efficiency;Background Rejection", 100, 0.0, 1.0);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.05);
    frame->Draw();

    // Draw ROC curves on top
    roc1->Draw("L SAME");
    roc2->Draw("L SAME");
    roc3->Draw("L SAME");
    ///*
    TMarker *pt = new TMarker(0.1951, 0.9962, 20);  // (x, y, marker style)
    pt->SetMarkerColor(kGreen+2);
    pt->SetMarkerSize(1);
    pt->Draw("SAME");
//*/
    // Add legend
    TLegend *legend = new TLegend(0.15, 0.15, 0.45, 0.3);
    legend->AddEntry(roc1, "BDT from Setup 1", "l");
    legend->AddEntry(roc2, "BDT from Setup 2", "l");
    legend->AddEntry(roc3, "BDT from Setup 3", "l");
    legend->Draw();

    gPad->Modified();
    gPad->Update();

    c1->SaveAs("ROC_Comparison_5.png");
}

