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

void FilterUtilities::Filter(double *waveform, double *filteredWaveform) {
    unsigned int GSFilterIdx = 0;
    unsigned int FFTFilterIdx = 0;
    intermediateInput = (double *) malloc(kNumberOfSamplingPoints * sizeof(double));
    intermediateOutput = (double *) malloc(kNumberOfSamplingPoints * sizeof(double));
    std::copy(waveform, waveform + kNumberOfSamplingPoints, intermediateInput);
    for (int i = 0; i < FilterSequence.size(); i++) {
        switch (FilterSequence.at(i)) {
            case kGS: {
                GSFilterCollection.at(GSFilterIdx)->Filter(intermediateInput, intermediateOutput,
                                                           kNumberOfSamplingPoints);
                std::copy(intermediateOutput, intermediateOutput + kNumberOfSamplingPoints, intermediateInput);
                GSFilterIdx++;
                break;
            }
            case kFFT: {
                FFTCutoffCollection.at(FFTFilterIdx)->Filter(intermediateInput, intermediateOutput);
                std::copy(intermediateOutput, intermediateOutput + kNumberOfSamplingPoints, intermediateInput);
                FFTFilterIdx++;
                break;
            }
            default: {
                std::cout << "Unknown filter types." << std::endl;
                break;
            }
        }
    }

    std::copy(intermediateOutput, intermediateOutput + kNumberOfSamplingPoints, filteredWaveform);
}

void FilterUtilities::AddFFTCutoffFilter(FreqRange cutRange) {
    FFTCutoff *fftCutoff = new FFTCutoff(kNumberOfSamplingPoints);
    fftCutoff->SetFreqCutoffRange(cutRange);
    FFTCutoffCollection.push_back(fftCutoff);
    FilterSequence.push_back(kFFT);
}

void FilterUtilities::SetNumberOfSamplingPoints(unsigned int numberOfSamplingPoints) {
    kNumberOfSamplingPoints = numberOfSamplingPoints;
}

void FilterUtilities::AddFFTCutoffFilterByValue(unsigned int freqCutValue) {
    FFTCutoff *fftCutoff = new FFTCutoff(kNumberOfSamplingPoints);
    fftCutoff->SetFreqCutoffValue(freqCutValue);
    FFTCutoffCollection.push_back(fftCutoff);
    FilterSequence.push_back(kFFT);
}






