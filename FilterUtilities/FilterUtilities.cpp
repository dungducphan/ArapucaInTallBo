//
// Created by Dung Phan on 10/3/17.
//

#include "FilterUtilities.h"

FilterUtilities::FilterUtilities() {
}

FilterUtilities::~FilterUtilities() {}

void FilterUtilities::AddGolaySavitzkyFilter(unsigned int NumberOfPoints, unsigned int Degree) {
    GSFilterCollection.push_back(new GolaySavitzkyCoeff(NumberOfPoints, Degree));
    FilterSequence.push_back(kGS);
}

void FilterUtilities::Filter(double *waveform, double *filteredWaveform, unsigned int sampleIdx) {
    unsigned int GSFilterIdx = 0;
    unsigned int FFTFilterIdx = 0;
    intermediateInput   = (double*)malloc(sampleIdx * sizeof(double));
    intermediateOutput  = (double*)malloc(sampleIdx * sizeof(double));
    std::copy(waveform, waveform + sampleIdx, intermediateInput);
    for (int i = 0; i < FilterSequence.size(); i++) {
        switch (FilterSequence.at(i)) {
            case kGS: {
                GSFilterCollection.at(GSFilterIdx)->Filter(intermediateInput, intermediateOutput, sampleIdx);
                std::copy(intermediateOutput, intermediateOutput + sampleIdx, intermediateInput);
                GSFilterIdx++;
                break;
            }
            case kFFT: {
                FFTFilterIdx++;
                break;
            }
            default: {
                std::cout << "Unknown filter types." << std::endl;
                break;
            }
        }
    }

    std::copy(intermediateOutput, intermediateOutput + sampleIdx, filteredWaveform);
}






