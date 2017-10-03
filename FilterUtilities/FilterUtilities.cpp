//
// Created by Dung Phan on 10/3/17.
//

#include "FilterUtilities.h"


FilterUtilities::FilterUtilities(unsigned int userNumberOfPoints, unsigned int userExtrapolationDegree) {
    kNumberOfPoints = userNumberOfPoints;
    kArrayMargin = (kNumberOfPoints - 1) / 2;
    kExtrapolationDegree = userExtrapolationDegree;
    GSFilter = new GolaySavitzkyCoeff(kNumberOfPoints, kExtrapolationDegree);

    kFuncSmoothingCoeffs = (double*)malloc(kNumberOfPoints * sizeof(double));
    GSFilter->GetFuncSmoothingCoeffs(kFuncSmoothingCoeffs);

    kFirstDerivativeSmoothingCoeffs = (double*)malloc(kNumberOfPoints * sizeof(double));
    GSFilter->GetFirstDerivativeSmoothingCoeffs(kFirstDerivativeSmoothingCoeffs);
}

FilterUtilities::~FilterUtilities() {}

void FilterUtilities::Filter(double *Waveform, double *FilteredWaveform, unsigned int NumberOfSamplingPoints) {
    for (int i = kArrayMargin; i < NumberOfSamplingPoints - kArrayMargin; i++) {
        FilteredWaveform[i] = ApplyFilteringWindow(Waveform, i);
    }

    for (int i = 0; i < kArrayMargin; i++) {
        FilteredWaveform[i] = FilteredWaveform[kArrayMargin];
    }

    for (int i = NumberOfSamplingPoints - kArrayMargin; i < NumberOfSamplingPoints; i++) {
        FilteredWaveform[i] = FilteredWaveform[NumberOfSamplingPoints - kArrayMargin - 1];
    }
}

double FilterUtilities::ApplyFilteringWindow(double *waveform, unsigned int sampleIdx) {
    unsigned int startingIdx = sampleIdx - kArrayMargin;
    double FilterValueAtSampleIdx = 0;
    for (int i = 0; i < kNumberOfPoints; i++) {
        FilterValueAtSampleIdx += waveform[startingIdx + i] * kFuncSmoothingCoeffs[i];
    }

    return FilterValueAtSampleIdx;
}



