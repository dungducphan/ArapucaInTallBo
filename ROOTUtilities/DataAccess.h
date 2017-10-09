#ifndef DataAccess_h
#define DataAccess_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>

class DataAccess {
public :
    TTree *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t fCurrent; //!current Tree number in a TChain
    TFile *f;

    // Declaration of leaf types
    Int_t adc[1650];

    // List of branches
    TBranch *b_adc;   //!

    DataAccess(TTree *tree = 0);

    virtual ~DataAccess();

    virtual Int_t Cut(Long64_t entry);

    virtual Int_t GetEntry(Long64_t entry);

    virtual Long64_t LoadTree(Long64_t entry);

    virtual void Init(TTree *tree);

    virtual void Loop();

    virtual Bool_t Notify();

    virtual void Show(Long64_t entry = -1);

    virtual Long64_t GetNEntries();

    virtual void AccessEntry(Long64_t ientry, double *accessPointer);
};

#endif