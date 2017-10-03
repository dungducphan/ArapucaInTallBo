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
    /*
     * This class performs Golay-Savitzky filter.
     * See more: https://en.wikipedia.org/wiki/Savitzky–Golay_filter.
     *
     * HOW TO USE:
     *
     *      - Create a GolaySavitzkyCoeff object.
     *      - Specify the number of points in a filtering window and the degree of filter polynomial.
     *      - To filter a waveform, call Filter(input, output, N).
     *      - The argument 'input' and 'output' are double array of N elements.
     *      - To print the smoothing matrices, call PrintMatrices().
     *      - To interface with a different class, call GetFuncSmoothingCoeffs().
     *
     * Example usage:
     *
     *      GolaySavitzkyCoeff* GSfilter = new GolaySavitzkyCoeff(5, 4);
     *      GSfilter->Filter(adc, filteredadc);
     *      GSfilter->PrintMatrices();
     *      GSfilter->GetFuncSmoothingCoeffs(outputCoeffs);
     */
    GolaySavitzkyCoeff(unsigned int userNumberOfPoints, unsigned int userExtrapolationDegree);
    virtual  ~GolaySavitzkyCoeff();

    /*
     * Printing smoothing matrices for debug purposes.
     */
    virtual void PrintMatrices();

    /*
     * Calculating the function smoothing coefficients and return them
     * as an array back to accessPointer.
     *
     * Use this to interface GolaySavitzkyCoeff with another class.
     */
    virtual void GetFuncSmoothingCoeffs(double* accessPointer);

    /*
     * Calculating the first derivative smoothing coefficients and
     * return them as an array back to accessPointer.
     *
     * Use this to interface GolaySavitzkyCoeff with another class.
     */
    virtual void GetFirstDerivativeSmoothingCoeffs(double* accessPointer);

    /*
     * Smoothing a waveform (array of values) using Golay-Savitzky
     * method.
     *
     *      - waveform          : input array of 'sampleIdx' double elements.
     *      - filteredWaveform  : input array of 'sampleIdx' double elements.
     *      - sampleIdx         : number of elements in a waveform.
     */
    virtual void Filter(double* waveform, double* filteredWaveform, unsigned int sampleIdx);

protected:

    /* -- protected set methods -- */
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