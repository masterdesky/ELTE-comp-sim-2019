#ifndef CPL_FFT_HPP
#define CPL_FFT_HPP

#include <complex>
#include <iostream>

namespace cpl {

class Vector;

class ComplexVector {
  public:
    ComplexVector(int dim = 1);
    ComplexVector(const ComplexVector& cv);
    ~ComplexVector() { delete [] v; }

    int dimension() const { return dim; }
    const std::complex<double> operator[](const int i) const { return v[i]; }
    std::complex<double>& operator[](const int i) { return v[i]; }
    ComplexVector& operator = (const ComplexVector& cv);

  private:
    int dim;
    std::complex<double> *v;
};

class FFT {
  public:
    FFT() { N = 0; f = 0; inverse = false; }
    void transform(ComplexVector& data);
    void inverseTransform(ComplexVector& data);
    Vector power(ComplexVector& data);

  private:
    int N;
    ComplexVector *f;
    bool inverse;

    void bitReverse();
    void DanielsonLanczos(int n);
};

class ComplexMatrix {
  public:

    ComplexMatrix(int rows=1, int cols=1, complex<double> d=0);
    ComplexMatrix(const ComplexMatrix& m);
    ~ComplexMatrix();

    int numRows() const { return rows; }
    int numCols() const { return cols; }
    const std::complex<double>* operator[](const int row) const;
    std::complex<double>* operator[](const int row);
    std::complex<double>& operator()(int row, int col);
    ComplexMatrix& operator=(const ComplexMatrix& m);

  private:
    int rows;
    int cols;
    std::complex<double> **m;
};

class FFT2 {
  public:
    FFT2() { }
    void transform(ComplexMatrix& data);
    void inverseTransform(ComplexMatrix& data);

  private:

};

}  /* end namespace cpl */

#endif /* CPL_FFT_HPP */
