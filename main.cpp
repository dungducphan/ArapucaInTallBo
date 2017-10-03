#include <iostream>

#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>

#include "ROOTUtilities/DataAccess.h"
#include "FilterUtilities/FilterUtilities.h"

int main() {
    DataAccess *myAnalysis_1 = new DataAccess();
    double *adc = (double *) malloc(1650 * sizeof(double));
    double* filteredadc = (double*)malloc(1650 * sizeof(double));
    double *idx = (double *) malloc(1650 * sizeof(double));
    myAnalysis_1->AccessEntry(40, adc);
    for (int i = 0; i < 1650; i++) {
        idx[i] = (double)i;
    }

    FilterUtilities* myFilter = new FilterUtilities();
    myFilter->SetNumberOfSamplingPoints(1650);
    myFilter->AddGolaySavitzkyFilter(5, 3);
    myFilter->AddFFTCutoffFilterByValue(800);
    myFilter->AddGolaySavitzkyFilter(25, 5);
    myFilter->Filter(adc, filteredadc);


    TGraph *adcGr = new TGraph(1650, idx, adc);
    TGraph* ftadcGr = new TGraph(1650, idx, filteredadc);
    ftadcGr->SetLineColor(kRed);
    ftadcGr->SetLineWidth(3);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(adcGr);
    mg->Add(ftadcGr);
    TCanvas* c = new TCanvas();
    mg->Draw("APL");
    c->SaveAs("testFilter.png");

    return 0;
}
