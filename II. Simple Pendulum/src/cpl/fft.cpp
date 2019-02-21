#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

#include "fft.hpp"
#include "vector.hpp"

namespace cpl {

ComplexVector::ComplexVector(int dim) {
    v = new std::complex<double> [this->dim = dim];
    for (int i = 0; i < dim; i++) v[i] = 0.0;
}

ComplexVector::ComplexVector(const ComplexVector& cv) {
    v = new std::complex<double> [dim = cv.dim];
    for (int i = 0; i < dim; i++) v[i] = cv.v[i];
}

ComplexVector& ComplexVector::operator = (const ComplexVector& cv) {
    if (this != &cv) {
        if (dim != cv.dim) {
            delete [] v;
            v = new std::complex<double> [dim = cv.dim];
        }
        for (int i = 0; i < dim; i++) v[i] = cv[i];
    }
    return *this;
}

// FFT implementation

void FFT::transform(ComplexVector& data) {
    N = data.dimension();
    f = &data;
    bitReverse();
    for (int n = 1; n < N; n *= 2)
        DanielsonLanczos(n);
    for (int i = 0; i < N; ++i)
        (*f)[i] /= std::sqrt(double(N));
}

void FFT::inverseTransform(ComplexVector& data) {
    inverse = true;
    transform(data);
    inverse = false;
}

void FFT::bitReverse() {
    int j = 1;
    for (int i = 1; i < N; ++i) {
        if (i < j) {
            std::complex<double> temp = (*f)[i-1];
            (*f)[i-1] = (*f)[j-1];
            (*f)[j-1] = temp;
        }
        int k = N / 2;
        while ( k < j ) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

void FFT::DanielsonLanczos(int n) {
    const double pi = 4 * atan(1.0);
    std::complex<double> W(0, pi / n);
    W = inverse ? std::exp(W) : std::exp(-W);
    std::complex<double> W_j(1, 0);
    for (int j = 0; j < n; ++j) {
        for (int i = j; i < N; i += 2 * n) {
            std::complex<double> temp = W_j * (*f)[n+i];
            (*f)[n+i] = (*f)[i] - temp;
            (*f)[i] += temp;
        }
        W_j *= W;
    }
}

Vector FFT::power(ComplexVector& data) {
    Vector P(1 + N / 2);
    P[0] = std::norm(data[0]) / double(N);
    for (int i = 1; i < N / 2; i++)
        P[i] = (std::norm(data[i]) + std::norm(data[N-i])) / double(N);
    P[N/2] = std::norm(data[N/2]) / double(N);
    return P;
}

static void error(const std::string str) {
    std::cerr << "cpl::ComplexMatrix error: " << str << std::endl;
    std::exit(EXIT_FAILURE);
}

static void error(const std::string str, int i) {
    std::ostringstream os;
    os << str << " " << i;
    error(os.str());
}

ComplexMatrix::ComplexMatrix(int rows, int cols, complex<double> d) {
    if (rows <= 0)
        error("bad number of rows", rows);
    if (cols <= 0)
        error("bad number of columns", cols);
    this->rows = rows;
    this->cols = cols;
    m = new complex<double>* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new complex<double> [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = d;
    }
}

ComplexMatrix::ComplexMatrix(const ComplexMatrix& mat) {
    rows = mat.rows;
    cols = mat.cols;

    m = new complex<double>* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new complex<double> [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = mat[i][j];
    }
}

ComplexMatrix::~ComplexMatrix() {
    for (int i = 0; i < rows; i++)
        delete [] m[i];
    delete [] m;
}

const complex<double>* ComplexMatrix::operator[](const int row) const {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}

complex<double>* ComplexMatrix::operator[](const int row) {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}

ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols) {
            for (int i = 0; i < rows; i++)
                delete [] m[i];
            delete [] m;
            rows = mat.rows;
            cols = mat.cols;
            m = new complex<double>* [rows];
            for (int i = 0; i < rows; i++)
                m[i] = new complex<double> [cols];
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] = mat[i][j];
    }
    return *this;
}

void FFT2::transform(ComplexMatrix& data) {
    FFT fft;
    int rows = data.numRows();
    int cols = data.numCols();
    ComplexVector r(cols), c(rows);
    // FFT rows of data
    for (int j = 0; j < rows; j++) {
        for (int k = 0; k < cols; k++)
            r[k] = data[j][k];
        fft.transform(r);
        for (int k = 0; k < cols; k++)
            data[j][k] = r[k];
    }
    // FFT columns of data
    for (int k = 0; k < cols; k++) {
        for (int j = 0; j < rows; j++)
            c[j] = data[j][k];
        fft.transform(c);
        for (int j = 0; j < rows; j++)
            data[j][k] = c[j];
    }
}

void FFT2::inverseTransform(ComplexMatrix& data) {
    FFT fft;
    int rows = data.numRows();
    int cols = data.numCols();
    ComplexVector r(cols), c(rows);
    // FFT rows of data
    for (int j = 0; j < rows; j++) {
        for (int k = 0; k < cols; k++)
            r[k] = data[j][k];
        fft.inverseTransform(r);
        for (int k = 0; k < cols; k++)
            data[j][k] = r[k];
    }
    // FFT columns of data
    for (int k = 0; k < cols; k++) {
        for (int j = 0; j < rows; j++)
            c[j] = data[j][k];
        fft.inverseTransform(c);
        for (int j = 0; j < rows; j++)
            data[j][k] = c[j];
    }
}

}  /* end namespace cpl */
