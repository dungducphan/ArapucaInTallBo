//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FFTCUTOFF_H
#define TALLBOANALYSIS_FFTCUTOFF_H

#include <fftw3.h>

class FFTCutoff {
public:
    FFTCutoff(unsigned int sampleIdx);
    ~FFTCutoff();

protected:
    unsigned int kNumberOfSamplingPoints;

private:
    fftw_complex *input, *filteredOutput;
    fftw_plan p;
};


#endif //TALLBOANALYSIS_FFTCUTOFF_H
