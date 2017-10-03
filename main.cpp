#include <iostream>
#include <armadillo>
#include <vector>

#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>

#include "ROOTUtilities/Analysis.h"
#include "FilterUtilities/GolaySavitzkyCoeff.h"
#include "FilterUtilities/FilterUtilities.h"

int main() {
    Analysis* myAnalysis_1 = new Analysis();
    Int_t* adc = (Int_t*)malloc(1650 * sizeof(Int_t));
    double* f_adc = (double*)malloc(1650 * sizeof(double));
    double* filteredadc = (double*)malloc(1650 * sizeof(double));
    double* dfilteredadc = (double*)malloc(1650 * sizeof(double));
    double* idx = (double*)malloc(1650 * sizeof(double));
    myAnalysis_1->AccessEntry(7, adc);
    for (int i = 0; i < 1650; i++) {
        idx[i] = (double)i;
        f_adc[i] = (double)adc[i];
    }

    FilterUtilities* myFilter = new FilterUtilities(17, 4);
    myFilter->Filter(f_adc, filteredadc, 1650);
    myFilter->Filter(filteredadc, dfilteredadc, 1650);
    myFilter->Filter(dfilteredadc, filteredadc, 1650);

    TGraph* adcGr = new TGraph(1650, idx, f_adc);
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