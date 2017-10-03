//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FFTCUTOFF_H
#define TALLBOANALYSIS_FFTCUTOFF_H

#include <iostream>
#include <fftw3.h>

enum FreqRange {kMedium, kHigh, kVeryHigh};

class FFTCutoff {
public:
    explicit FFTCutoff(unsigned int sampleIdx);
    ~FFTCutoff();

    virtual void Filter(double* waveform, double* filteredWaveform);
    virtual void SetFreqCutoffRange(FreqRange cutRange);
    virtual void SetFreqCutoffValue(unsigned int freqCutValue);

protected:
    virtual void ForwardFFT(double *waveform);
    virtual void BackwardFFT(double* filteredWaveform);
    virtual void CutOffHighFreqNoise();
    virtual void DestroyPlan();

private:
    fftw_complex *input;
    fftw_complex *filteredOutput;
    fftw_complex *intermediateSpectrum;
    fftw_plan p_forward;
    fftw_plan p_backward;

    FreqRange kCutRange;
    unsigned int kCutValue;
    bool kCutByRange;

    unsigned int kSampleToStartCutting;
    unsigned int kNumberOfSamplingPoints;
};


#endif //TALLBOANALYSIS_FFTCUTOFF_H
