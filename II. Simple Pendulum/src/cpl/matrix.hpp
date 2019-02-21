#ifndef CPL_MATRIX_HPP
#define CPL_MATRIX_HPP

namespace cpl {

class Vector;

class Matrix {
  public:

    Matrix(int rows=1, int cols=1, double d=0);

    Matrix(const Matrix& m);

    ~Matrix();

    int numRows() const { return rows; }

    int numCols() const { return cols; }

    const double* operator[](const int row) const;

    double* operator[](const int row);

    double& operator()(int row, int col);

    Matrix& operator=(const Matrix& m);

    Matrix& operator=(const double d);

    Matrix& operator += (const Matrix& m);

    Matrix& operator -= (const Matrix& m);

    Matrix& operator *= (double d);

    Matrix& operator /= (double d);

    Matrix transpose();

    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);

  private:
    int rows;
    int cols;
    double **m;
};

inline Matrix operator + (const Matrix& m) {
    return m;
}

extern Matrix operator - (const Matrix& m);

extern Matrix operator * (const Matrix& m, double d);

extern Matrix operator * (double d, const Matrix& m);

extern Matrix operator / (const Matrix& m, double d);

extern Matrix operator + (const Matrix& m1, const Matrix& m2);

extern Matrix operator - (const Matrix& m1, const Matrix& m2);

extern Matrix operator * (const Matrix& m1, const Matrix& m2);

extern void solveGaussJordan(Matrix& A, Matrix& B);

extern void solveLUDecompose(Matrix& A, Matrix& B);

extern void reduceHouseholder(Matrix& A, Vector& d, Vector& e);

extern void solveTQLI(Vector& d, Vector&e, Matrix& z);

extern void sortEigen(Vector& d, Matrix& z);

extern Vector solveEigenSymmetric(Matrix& A);

extern void singularValueDecompose(Matrix& a, Vector& w, Matrix& v);

extern void minimizeBFGS(Vector& p, const double gtol, int& iter, double& fret,
            double (*func)(Vector&), void (*dfunc)(Vector&, Vector&));

} /* end namespace cpl */

#endif /* CPL_MATRIX_HPP */
