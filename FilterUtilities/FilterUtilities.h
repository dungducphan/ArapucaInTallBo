//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FILTERUTILITIES_H
#define TALLBOANALYSIS_FILTERUTILITIES_H

#include "GolaySavitzkyCoeff.h"

enum FilterTypes {kGS, kFFT};

class FilterUtilities {
public:
    FilterUtilities();
    virtual ~FilterUtilities();

    virtual void AddGolaySavitzkyFilter(unsigned int NumberOfPoints, unsigned int Degree);
    virtual void Filter(double* waveform, double* filteredWaveform, unsigned int sampleIdx);

private:
    double* intermediateInput;
    double* intermediateOutput;
    std::vector<GolaySavitzkyCoeff*> GSFilterCollection;
    std::vector<FilterTypes> FilterSequence;
};

#endif //TALLBOANALYSIS_FILTERUTILITIES_H
