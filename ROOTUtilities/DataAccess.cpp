#include "DataAccess.h"

DataAccess::DataAccess(TTree *tree) : fChain(0) {
    if (tree == 0) {
        TFile *f = (TFile *) gROOT->GetListOfFiles()->FindObject(
                "../FreeRun-2017-09-22-16-32-Threshold_15-SubRun-ch0.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("../FreeRun-2017-09-22-16-32-Threshold_15-SubRun-ch0.root");
        }
        f->GetObject("waveform", tree);

    }
    Init(tree);
}

DataAccess::~DataAccess() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t DataAccess::GetEntry(Long64_t entry) {
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t DataAccess::LoadTree(Long64_t entry) {
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void DataAccess::Init(TTree *tree) {
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("adc", adc, &b_adc);
    Notify();
}

Bool_t DataAccess::Notify() {
    return kTRUE;
}

void DataAccess::Show(Long64_t entry) {
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t DataAccess::Cut(Long64_t entry) {
// returns  1 if entry is accepted.
// returns -1 otherwise.
    return 1;
}

void DataAccess::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
    }
}

Long64_t DataAccess::GetNEntries() {
    return fChain->GetEntriesFast();
}

void DataAccess::AccessEntry(Long64_t ientry, double *accessPointer) {
    Long64_t ientryChain = LoadTree(ientry);
    Long64_t nb = fChain->GetEntry(ientry);
    for (int i = 0; i < 1650; i++) {
        accessPointer[i] = (double) adc[i];
    }
}


