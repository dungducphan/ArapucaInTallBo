#ifndef GSFilter_h
#define GSFilter_h

#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TMath.h>

class GolaySavitzkyCoeff {
public:
    GolaySavitzkyCoeff(unsigned int userNumberOfPoints, unsigned int userExtrapolationDegree);
    virtual  ~GolaySavitzkyCoeff();
    virtual void PrintMatrices();

    // interfacing methods
    virtual void GetFuncSmoothingCoeffs(double* accessPointer);
    virtual void GetFirstDerivativeSmoothingCoeffs(double* accessPointer);
    virtual void Filter(double* waveform, double* filteredWaveform, unsigned int sampleIdx);

protected:

    // protected set
    virtual void SetNumberOfPoints(unsigned int);
    virtual void SetExtrapolationDegree(unsigned int);

    // protected implementation
    virtual double ApplyFilteringWindow(double* waveform, unsigned int sampleIdx);

    // protected data members
    unsigned int kNumberOfPoints;
    unsigned int kArrayMargin;
    unsigned int kExtrapolationDegree;

    // protected interfacings data (for derived classes)
    double*   kFuncSmoothingCoeffs;
    double*   kFirstDerivativeSmoothingCoeffs;

private:

    // private implementations
    virtual void CalculateJMatrix();
    virtual void CalculateAEvenMatrix();
    virtual void CalculateAOddMatrix();
    virtual void CalculateCMatrix();
    virtual void CalculateFuncSmoothingCoeffs();
    virtual void CalculateFirstDerivativeSmoothingCoeffs();

    // private data members
    double*   JArray;
    double*   JJTransposeArray;
    TMatrixD* JMatrix;
    TMatrixD* JTransposeMatrix;
    TMatrixD* JJTransposeMatrix;
    TMatrixD* CMatrix;

    TMatrixD* AEvenMatrix;
    TMatrixD* AOddMatrix;
    TMatrixD* JEvenMatrix;
    TMatrixD* JJTransposeEvenMatrix;
    TMatrixD* JOddMatrix;
    TMatrixD* JJTransposeOddMatrix;
};

#endif