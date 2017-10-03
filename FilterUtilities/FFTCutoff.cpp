//
// Created by Dung Phan on 10/3/17.
//

#include "FFTCutoff.h"

FFTCutoff::FFTCutoff(unsigned int sampleIdx) {
    kNumberOfSamplingPoints = sampleIdx;
    kCutByRange = true;
    input                   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * kNumberOfSamplingPoints);
    filteredOutput          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * kNumberOfSamplingPoints);
    intermediateSpectrum    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * kNumberOfSamplingPoints);
}

FFTCutoff::~FFTCutoff() {
}

void FFTCutoff::Filter(double *waveform, double *filteredWaveform) {
    ForwardFFT(waveform);
    CutOffHighFreqNoise();
    BackwardFFT(filteredWaveform);
    DestroyPlan();
}

void FFTCutoff::ForwardFFT(double *waveform) {
    for (int i = 0; i < kNumberOfSamplingPoints; i++) {
        input[i][0] = waveform[i];
        input[i][1] = 0;
    }
    p_forward = fftw_plan_dft_1d(kNumberOfSamplingPoints, input, intermediateSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);
}

void FFTCutoff::BackwardFFT(double *filteredWaveform) {
    p_backward = fftw_plan_dft_1d(kNumberOfSamplingPoints, intermediateSpectrum, filteredOutput, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_backward);
    for (int i = 0; i < kNumberOfSamplingPoints; i++) {
        // FFT scaling problem
        // Need to divide by 'kNumberOfSamplingPoints'
        // http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs
        filteredWaveform[i] = filteredOutput[i][0]/kNumberOfSamplingPoints;
    }
}

void FFTCutoff::SetFreqCutoffRange(FreqRange cutRange) {
    kCutRange = cutRange;
}

void FFTCutoff::CutOffHighFreqNoise() {
    if (kCutByRange) {
        switch(kCutRange) {
            case kMedium: {
                kSampleToStartCutting = (kNumberOfSamplingPoints * 1) / 4;
                break;
            }
            case kHigh: {
                kSampleToStartCutting = (kNumberOfSamplingPoints * 2) / 4;
                break;
            }
            case kVeryHigh: {
                kSampleToStartCutting = (kNumberOfSamplingPoints * 3) / 4;
                break;
            }
            default: {
                std::cout << "Unknown frequency cut." << std::endl;
            }
        }

        for (int i = kSampleToStartCutting; i < kNumberOfSamplingPoints; i++) {
            intermediateSpectrum[i][0] = 0;
        }
    } else {
        for (int i = kCutValue; i < kNumberOfSamplingPoints; i++) {
            intermediateSpectrum[i][0] = 0;
        }
    }
}

void FFTCutoff::DestroyPlan() {
    fftw_destroy_plan(p_backward);
    fftw_destroy_plan(p_forward);
    fftw_free(input);
    fftw_free(intermediateSpectrum);
    fftw_free(filteredOutput);
}

void FFTCutoff::SetFreqCutoffValue(unsigned int freqCutValue) {
    if (freqCutValue > kNumberOfSamplingPoints) {
        std::cout << "Freq cut value cannot exceed the number of sampling points." << std::endl;
    } else {
        kCutValue = freqCutValue;
        kCutByRange = false;
    }
}
