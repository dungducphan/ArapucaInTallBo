//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FILTERUTILITIES_H
#define TALLBOANALYSIS_FILTERUTILITIES_H

#include "GolaySavitzkyCoeff.h"
#include "FFTCutoff.h"

enum FilterTypes {kGS, kFFT};

class FilterUtilities {
public:
    FilterUtilities();
    virtual ~FilterUtilities();

    virtual void AddGolaySavitzkyFilter(unsigned int NumberOfPoints, unsigned int Degree);

    virtual void AddFFTCutoffFilter(FreqRange cutRange);

    virtual void AddFFTCutoffFilterByValue(unsigned int freqCutValue);

    virtual void Filter(double *waveform, double *filteredWaveform);

    virtual void SetNumberOfSamplingPoints(unsigned int numberOfSamplingPoints);

private:
    double* intermediateInput;
    double* intermediateOutput;
    std::vector<GolaySavitzkyCoeff*> GSFilterCollection;
    std::vector<FFTCutoff *> FFTCutoffCollection;
    std::vector<FilterTypes> FilterSequence;

    unsigned int kNumberOfSamplingPoints;
};

#endif //TALLBOANALYSIS_FILTERUTILITIES_H
