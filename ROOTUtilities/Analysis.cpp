#define Analysis_cxx

#include "Analysis.h"

Analysis::Analysis(TTree *tree) : fChain(0) {
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

Analysis::~Analysis() {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry) {
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t Analysis::LoadTree(Long64_t entry) {
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void Analysis::Init(TTree *tree) {
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("adc", adc, &b_adc);
    Notify();
}

Bool_t Analysis::Notify() {
    return kTRUE;
}

void Analysis::Show(Long64_t entry) {
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t Analysis::Cut(Long64_t entry) {
// returns  1 if entry is accepted.
// returns -1 otherwise.
    return 1;
}

void Analysis::Loop() {
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

Long64_t Analysis::GetNEntries() {
    return fChain->GetEntriesFast();
}

void Analysis::AccessEntry(Long64_t ientry, Int_t *accessPointer) {
    Long64_t ientryChain = LoadTree(ientry);
    Long64_t nb = fChain->GetEntry(ientry);
    std::copy(adc, adc + 1650, accessPointer);
}


