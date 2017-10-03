//
// Created by Dung Phan on 10/3/17.
//

#ifndef TALLBOANALYSIS_FILTERUTILITIES_H
#define TALLBOANALYSIS_FILTERUTILITIES_H

#include "GolaySavitzkyCoeff.h"
#include "FFTCutoff.h"

/*
 * This enum specifies the type of filters user want to perform
 *      - kGS   : Golay-Savitzky filter.
 *      - kFFT  : FFT filter.
 */
enum FilterTypes {kGS, kFFT};

class FilterUtilities {
public:
    /*
     * This class is a wrapper that allows user to choose
     * to perform different types of filters on the data,
     * as well as repeat a specific filter many times.
     *
     * It provides methods to create, configure and add
     * wanted filters into the filter sequence to be performed.
     *
     * HOW TO USE:
     *      - Create a FilterUtilities object.
     *      - Set number of sampling points N of the data.
     *      - Add wanted filters method with their parameters. See their classes for details.
     *      - Call Filter(input, output) function to perform filter sequence.
     *      - Notice that input and output have to be double array with N elements.
     *
     * Example usage:
     *
     *      FilterUtilities* myFilter = new FilterUtilities();
     *      myFilter->SetNumberOfSamplingPoints(1650);
     *      myFilter->AddGolaySavitzkyFilter(5, 3);
     *      myFilter->AddFFTCutoffFilterByValue(800);
     *      myFilter->AddGolaySavitzkyFilter(25, 5);
     *      myFilter->Filter(adc, filteredadc);
     */
    FilterUtilities();
    virtual ~FilterUtilities();

    /*
     * Set the number of sampling points in a waveform.
     *
     * For e.g.: If waveform contains of 1024 ADC values, then we call:
     *
     *      FilterUtilities* myFilter = new FilterUtilities();
     *      myFilter->SetNumberOfSamplingPoints(1024);
     */
    virtual void SetNumberOfSamplingPoints(unsigned int numberOfSamplingPoints);

    /*
     * IMPORTANT NOTICE:
     * User needs to call the Add...Filter() methods in the order in which
     * they want the filters are perform.
     *
     * For e.g.: To perform a sequence of Golay-Savitzky filter, then
     *           FFT filter and then another Golay-Savitzky, we call:
     *
     *      FilterUtilities* myFilter = new FilterUtilities();
     *      myFilter->SetNumberOfSamplingPoints(1650);
     *      myFilter->AddGolaySavitzkyFilter(5, 3);
     *      myFilter->AddFFTCutoffFilterByValue(800);
     *      myFilter->AddGolaySavitzkyFilter(25, 5);
     *
     * The order of filters to be performed is then saved in FilterSequence.
     */

    /*
     * Adding GolaySavitzkyFilter object
     * Need to specify number of points in a filter window and the degree of the polynomial filter
     */
    virtual void AddGolaySavitzkyFilter(unsigned int NumberOfPoints, unsigned int Degree);

    /*
     * Adding GolaySavitzkyFilter object
     * Need to specify the cutoff frequency
     * Cutoff range is pre-defined in the FFTCutoff class
     *      - kMedium   : cutoff starts at 1/4 of the waveform.
     *      - kHigh     : cutoff starts at 1/2 of the waveform.
     *      - kVeryHigh : cutoff starts at 3/4 of the waveform.
     */
    virtual void AddFFTCutoffFilter(FreqRange cutRange);

    /*
     * Adding GolaySavitzkyFilter object
     * Need to specify the cutoff frequency
     * Cutoff freq is an index number between 0 and kNumberOfSamplingPoints
     */
    virtual void AddFFTCutoffFilterByValue(unsigned int freqCutValue);

    /*
     * Perform filtering on data
     *      - waveform          : input waveform.
     *      - filteredWaveform  : filtered output waveform.
     *
     * Filter() calls different filtering methods in an user-defined order specified
     * by FilterSequence.
     */
    virtual void Filter(double *waveform, double *filteredWaveform);

private:
    double* intermediateInput;
    double* intermediateOutput;
    std::vector<GolaySavitzkyCoeff*> GSFilterCollection;
    std::vector<FFTCutoff *> FFTCutoffCollection;
    std::vector<FilterTypes> FilterSequence;

    /*
     * Number of sampling points of a waveform.
     */
    unsigned int kNumberOfSamplingPoints;
};

#endif //TALLBOANALYSIS_FILTERUTILITIES_H
