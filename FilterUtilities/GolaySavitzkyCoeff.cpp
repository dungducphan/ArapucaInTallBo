#include "GolaySavitzkyCoeff.h"


GolaySavitzkyCoeff::GolaySavitzkyCoeff(unsigned int userNumberOfPoints, unsigned int userExtrapolationDegree) {
    SetExtrapolationDegree(userExtrapolationDegree);
    SetNumberOfPoints(userNumberOfPoints);
    kArrayMargin = (kNumberOfPoints - 1) / 2;
    CalculateJMatrix();
    CalculateAEvenMatrix();
    CalculateAOddMatrix();
    CalculateCMatrix();
    CalculateFuncSmoothingCoeffs();
    CalculateFirstDerivativeSmoothingCoeffs();
}

GolaySavitzkyCoeff::~GolaySavitzkyCoeff() {
}

void GolaySavitzkyCoeff::PrintMatrices() {

    std::cout << "\t The J matrix: " << std::endl;
    JMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The J-transposed matrix: " << std::endl;
    JTransposeMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The product of J and J-transposed matrices: " << std::endl;
    JJTransposeMatrix->Print();
    std::cout << std::endl;

    std::cout << "The C matrix: " << std::endl;
    CMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The J-even matrices: " << std::endl;
    JEvenMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The JJTranspose-even matrices: " << std::endl;
    JJTransposeEvenMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The A-even matrices: " << std::endl;
    AEvenMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The J-odd matrices: " << std::endl;
    JOddMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The JJTranspose-odd matrices: " << std::endl;
    JJTransposeOddMatrix->Print();
    std::cout << std::endl;

    std::cout << "\t The A-odd matrices: " << std::endl;
    AOddMatrix->Print();
    std::cout << std::endl;

    std::cout << "The function smoothing coefficents using Golay-Savitzky filter are: " << std::endl;
    for (unsigned int i = 0; i < kNumberOfPoints; i++) {
        std::cout << "A[" << i << "] = " << std::setw(12) << *(kFuncSmoothingCoeffs + i) << "." << std::endl;
    }
    std::cout << std::endl;

    std::cout << "The first derivative smoothing coefficents using Golay-Savitzky filter are: " << std::endl;
    for (unsigned int i = 0; i < kNumberOfPoints; i++) {
        std::cout << "B[" << i << "] = " << std::setw(12) << *(kFirstDerivativeSmoothingCoeffs + i) << "." << std::endl;
    }
    std::cout << std::endl;

    return;
}

void GolaySavitzkyCoeff::SetNumberOfPoints(unsigned int userNumberOfPoints) {
    if (userNumberOfPoints % 2 == 1) {
        kNumberOfPoints = userNumberOfPoints;
    } else {
        std::cout << "\t WARNING: Number of interpolation points has to be an odd integer." << std::endl;
        std::cout << "\t WARNING: Filter algorithm automatically sets number of points to " << userNumberOfPoints + 1 << "." << std::endl;
        kNumberOfPoints = userNumberOfPoints + 1;
    }

    return;
}

void GolaySavitzkyCoeff::SetExtrapolationDegree(unsigned int userExtrapolationDegree) {
    kExtrapolationDegree = userExtrapolationDegree;
    return;
}

void GolaySavitzkyCoeff::CalculateJMatrix() {
    JArray = (double*) malloc(kNumberOfPoints * kExtrapolationDegree * sizeof(double));
    for (int i = 0; i < kNumberOfPoints; i++) {
        for (int j = 0; j < kExtrapolationDegree; j++) {
            if ((- ((kNumberOfPoints - 1) / 2) + i == 0) && (j == 0)) {
                *(JArray + kExtrapolationDegree * i + j) = 1;
            } else {
                *(JArray + kExtrapolationDegree * i + j) = TMath::Power(- (((double) kNumberOfPoints - 1) / 2) + (double) i, (double) j);
            }
        }
    }

    JMatrix = new TMatrixD(kNumberOfPoints, kExtrapolationDegree, JArray);

    JTransposeMatrix = new TMatrixD(kExtrapolationDegree, kNumberOfPoints);
    JTransposeMatrix->Transpose(*JMatrix);

    JJTransposeMatrix = new TMatrixD(kExtrapolationDegree, kExtrapolationDegree);
    JJTransposeMatrix->Mult(*JTransposeMatrix, *JMatrix);
    JJTransposeArray = JJTransposeMatrix->GetMatrixArray();

    return;
}

void GolaySavitzkyCoeff::CalculateAEvenMatrix() {
    unsigned int NumberOfColumns;
    if (kExtrapolationDegree % 2 == 0) {
        NumberOfColumns = (int)(kExtrapolationDegree / 2);
    } else {
        NumberOfColumns = (int)(kExtrapolationDegree / 2) + 1;
    }
    unsigned int NumberOfRows    = kNumberOfPoints;

    double* JEvenArray = (double*) malloc(NumberOfRows * NumberOfColumns * sizeof(double));
    for (unsigned int i = 0; i < NumberOfRows; i++) {
        for (unsigned int j = 0; j < kExtrapolationDegree; j++) {
            if (j % 2 == 0) {
                *(JEvenArray + i*NumberOfColumns + (int)(j/2)) = *(JArray + i*kExtrapolationDegree + j);
            }
        }
    }
    JEvenMatrix = new TMatrixD(NumberOfRows, NumberOfColumns, JEvenArray);

    double* JJTransposeEvenArray = (double*) malloc(NumberOfColumns * NumberOfColumns * sizeof(double));
    for (unsigned int i = 0; i < kExtrapolationDegree; i++) {
        for (unsigned int j = 0; j < kExtrapolationDegree; j++) {
            if (j % 2 == 0 && i % 2 == 0) {
                *(JJTransposeEvenArray + (int)(i/2)*NumberOfColumns + (int)(j/2)) = *(JJTransposeArray + i*kExtrapolationDegree + j);
            }
        }
    }
    JJTransposeEvenMatrix = new TMatrixD(NumberOfColumns, NumberOfColumns, JJTransposeEvenArray);

    /* Calculate AEven */
    // Inverse JJTransposeEven matrix
    double Determinant;
    JJTransposeEvenMatrix->InvertFast(&Determinant);
    // Transposer JEven matrix
    TMatrixD* JTransposeEvenMatrix = new TMatrixD(NumberOfColumns, NumberOfRows);
    JTransposeEvenMatrix->Transpose(*JEvenMatrix);
    // Important notice
    // - NumberOfRows happens to be the number of columns of the AEven matrix;
    // - NumberOfColumns happens to be the number of rows of the AEven matrix;
    // - AEven matrix is the product of Inversed-JJTransposeEven and JTransposeEven
    AEvenMatrix = new TMatrixD(NumberOfColumns, NumberOfRows);
    AEvenMatrix->Mult(*JJTransposeEvenMatrix, *JTransposeEvenMatrix);

    return;
}

void GolaySavitzkyCoeff::CalculateAOddMatrix() {
    unsigned int NumberOfColumns;
    if (kExtrapolationDegree % 2 == 0) {
        NumberOfColumns = (int)(kExtrapolationDegree / 2);
    } else {
        NumberOfColumns = (int)(kExtrapolationDegree / 2);
    }
    unsigned int NumberOfRows    = kNumberOfPoints;

    double* JOddArray = (double*) malloc(NumberOfRows * NumberOfColumns * sizeof(double));
    for (unsigned int i = 0; i < NumberOfRows; i++) {
        for (unsigned int j = 0; j < kExtrapolationDegree; j++) {
            if (j % 2 == 1) {
                *(JOddArray + i*NumberOfColumns + (int)(j/2)) = *(JArray + i*kExtrapolationDegree + j);
            }
        }
    }
    JOddMatrix = new TMatrixD(NumberOfRows, NumberOfColumns, JOddArray);

    double* JJTransposeOddArray = (double*) malloc(NumberOfColumns * NumberOfColumns * sizeof(double));
    for (unsigned int i = 0; i < kExtrapolationDegree; i++) {
        for (unsigned int j = 0; j < kExtrapolationDegree; j++) {
            if (j % 2 == 1 && i % 2 == 1) {
                *(JJTransposeOddArray + (int)(i/2)*NumberOfColumns + (int)(j/2)) = *(JJTransposeArray + i*kExtrapolationDegree + j);
            }
        }
    }
    JJTransposeOddMatrix = new TMatrixD(NumberOfColumns, NumberOfColumns, JJTransposeOddArray);

    /* Calculate AOdd */
    // Inverse JJTransposeOdd matrix
    double Determinant;
    JJTransposeOddMatrix->InvertFast(&Determinant);
    // Transposer JOdd matrix
    TMatrixD* JTransposeOddMatrix = new TMatrixD(NumberOfColumns, NumberOfRows);
    JTransposeOddMatrix->Transpose(*JOddMatrix);
    // Important notice
    // - NumberOfRows happens to be the number of columns of the AOdd matrix;
    // - NumberOfColumns happens to be the number of rows of the AOdd matrix;
    // - AOdd matrix is the product of Inversed-JJTransposeOdd and JTransposeOdd
    AOddMatrix = new TMatrixD(NumberOfColumns, NumberOfRows);
    AOddMatrix->Mult(*JJTransposeOddMatrix, *JTransposeOddMatrix);

    return;
}

void GolaySavitzkyCoeff::CalculateCMatrix() {
    double Determinant;
    JJTransposeMatrix->InvertFast(&Determinant);

    CMatrix = new TMatrixD(kExtrapolationDegree, kNumberOfPoints);
    CMatrix->Mult(*JJTransposeMatrix, *JTransposeMatrix);

    return;
}

void GolaySavitzkyCoeff::CalculateFuncSmoothingCoeffs() {
    kFuncSmoothingCoeffs = (double*) malloc(kNumberOfPoints * sizeof(double));
    AEvenMatrix->ExtractRow(0, 0, kFuncSmoothingCoeffs);

    return;
}

void GolaySavitzkyCoeff::CalculateFirstDerivativeSmoothingCoeffs() {
    kFirstDerivativeSmoothingCoeffs = (double*) malloc(kNumberOfPoints * sizeof(double));
    AOddMatrix->ExtractRow(0, 0, kFirstDerivativeSmoothingCoeffs);

    return;
}

void GolaySavitzkyCoeff::GetFuncSmoothingCoeffs(double *accessPointer) {
    std::copy(kFuncSmoothingCoeffs, kFuncSmoothingCoeffs + kNumberOfPoints, accessPointer);
}

void GolaySavitzkyCoeff::GetFirstDerivativeSmoothingCoeffs(double *accessPointer) {
    std::copy(kFirstDerivativeSmoothingCoeffs, kFirstDerivativeSmoothingCoeffs + kNumberOfPoints, accessPointer);
}

void GolaySavitzkyCoeff::Filter(double *Waveform, double *FilteredWaveform, unsigned int NumberOfSamplingPoints) {
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

double GolaySavitzkyCoeff::ApplyFilteringWindow(double *waveform, unsigned int sampleIdx) {
    unsigned int startingIdx = sampleIdx - kArrayMargin;
    double FilterValueAtSampleIdx = 0;
    for (int i = 0; i < kNumberOfPoints; i++) {
        FilterValueAtSampleIdx += waveform[startingIdx + i] * kFuncSmoothingCoeffs[i];
    }

    return FilterValueAtSampleIdx;
}

/*
void GolaySavitzkyCoeff::FirstDerivative(double* Waveform, double* FilteredFirstDerivative, double VariableStep) {
    int HalfNumberOfPoints = (int)(kNumberOfPoints / 2);
    for (unsigned int i = HalfNumberOfPoints; i < 1024 - HalfNumberOfPoints; i++) {
        *(FilteredFirstDerivative + i) = 0;
        for (int j = -HalfNumberOfPoints; j < HalfNumberOfPoints + 1; j++) {
            *(FilteredFirstDerivative + i) += *(Waveform + i + j) * *(kFirstDerivativeSmoothingCoeffs + j + HalfNumberOfPoints);
        }
        *(FilteredFirstDerivative + i) = *(FilteredFirstDerivative + i) / VariableStep;
    }

    for (unsigned int i = 0; i < HalfNumberOfPoints; i++) {
        *(FilteredFirstDerivative + i) = *(FilteredFirstDerivative + HalfNumberOfPoints);
        *(FilteredFirstDerivative + 1023 - i) = *(FilteredFirstDerivative + 1023 - HalfNumberOfPoints);
    }

    return;
}

*/