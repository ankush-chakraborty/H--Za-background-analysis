#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCut.h>
#include <TMVA/Factory.h>
#include <TMVA/Tools.h>
#include <TMVA/DataLoader.h>
#include <iostream>

void TMVAClassification_BDT3() {
    TMVA::Tools::Instance();
    
    TFile *outputFile = TFile::Open("TMVAOutput_3.root", "RECREATE");
    
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
    
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
    
    TFile *input = TFile::Open("tmva_input1.root");
    if (!input || input->IsZombie()) {
        std::cout << "Error: Cannot open input file tmva_input1.root" << std::endl;
        return;
    }

    TTree *tree = (TTree*)input->Get("tmva_tree");
    if (!tree) {
        std::cout << "Error: Cannot find tmva_tree in input file" << std::endl;
        input->Close();
        return;
    }
    
    dataloader->AddSignalTree(tree, 1.0);
    dataloader->AddBackgroundTree(tree, 1.0);
    
    dataloader->SetSignalWeightExpression("weight");
    dataloader->SetBackgroundWeightExpression("weight");
    
    dataloader->AddVariable("deta_jj", 'F');
    dataloader->AddVariable("dphi_jj", 'F');
    dataloader->AddVariable("dphi_jj_zg", 'F');
    dataloader->AddVariable("dR_j_g", 'F');
    dataloader->AddVariable("pT_j1", 'F');
    dataloader->AddVariable("pT_j2", 'F');
    dataloader->AddVariable("pTt", 'F');
    dataloader->AddVariable("pT_balance", 'F');
    //dataloader->AddVariable("zepp_g", 'F');
    dataloader->AddVariable("zepp_h", 'F');
    //dataloader->AddVariable("mjj", 'F');
    //dataloader->AddVariable("njets", 'I');
    dataloader->AddVariable("njets2", 'F');
    //dataloader->AddVariable("l1_pt", 'F');
    //dataloader->AddVariable("l2_pt", 'F');
    //dataloader->AddVariable("pho_pt", 'F');
    //dataloader->AddVariable("pho_eta", 'F');
    //dataloader->AddVariable("z_pt", 'F');
    //dataloader->AddVariable("l_eta", 'F');
    //dataloader->AddVariable("var1", 'F');
    //dataloader->AddVariable("var2", 'F');
    //dataloader->AddVariable("var3", 'F');
    
    TCut sigCut = "classID == 1";
    TCut bkgCut = "classID == 0";

    dataloader->PrepareTrainingAndTestTree(sigCut, bkgCut,
        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:!V");
    
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
        "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:"
        "AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");
    
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();
    input->Close();

    delete factory;
    delete dataloader;

    std::cout << "TMVA training complete! Results saved in TMVAOutput.root" << std::endl;
}

