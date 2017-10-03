//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FILTERUTILITIES_H
#define TALLBOANALYSIS_FILTERUTILITIES_H

#include "GolaySavitzkyCoeff.h"

class FilterUtilities {
public:
    FilterUtilities(unsigned int userNumberOfPoints, unsigned int userExtrapolationDegree);
    virtual ~FilterUtilities();

    virtual void Filter(double* waveform, double* filteredWaveform, unsigned int sampleIdx);
    // virtual void FirstDerivative(double* Waveform, double* FilteredWaveform, double VariableStep = 1.);

private:
    GolaySavitzkyCoeff* GSFilter;

    unsigned int kNumberOfPoints;
    unsigned int kArrayMargin;
    unsigned int kExtrapolationDegree;

    double*   kFuncSmoothingCoeffs;
    double*   kFirstDerivativeSmoothingCoeffs;

    // Private methods
    virtual double ApplyFilteringWindow(double* waveform, unsigned int sampleIdx);
};

#endif //TALLBOANALYSIS_FILTERUTILITIES_H
