#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>

#include "matrix.hpp"
#include "vector.hpp"

namespace cpl {

static void error(const std::string str) {
    std::cerr << "cpl::Matrix error: " << str << std::endl;
    std::exit(EXIT_FAILURE);
}

static void error(const std::string str, int i) {
    std::ostringstream os;
    os << str << " " << i;
    error(os.str());
}

Matrix::Matrix(int rows, int cols, double d) {
    if (rows <= 0)
        error("bad number of rows", rows);
    if (cols <= 0)
        error("bad number of columns", cols);
    this->rows = rows;
    this->cols = cols;
    m = new double* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new double [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = d;
    }
}

Matrix::Matrix(const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;

    m = new double* [rows];
    for (int i = 0; i < rows; i++) {
        m[i] = new double [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = mat[i][j];
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < rows; i++)
        delete [] m[i];
    delete [] m;
}

const double* Matrix::operator[](const int row) const {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}

double* Matrix::operator[](const int row) {
    if (row < 0 || row >= rows)
        error("bad row index", row);
    return m[row];
}

Matrix& Matrix::operator=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols) {
            for (int i = 0; i < rows; i++)
                delete [] m[i];
            delete [] m;
            rows = mat.rows;
            cols = mat.cols;
            m = new double* [rows];
            for (int i = 0; i < rows; i++)
                m[i] = new double [cols];
        }
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] = mat[i][j];
    }
    return *this;
}

Matrix& Matrix::operator=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] = d;
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols)
            error("matrix dimension mismatch");
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] += mat[i][j];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& mat) {
    if (this != &mat) {
        if (rows != mat.rows || cols != mat.cols)
            error("matrix dimension mismatch");
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] -= mat[i][j];
    }
    return *this;
}

Matrix& Matrix::operator*=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] *= d;
    return *this;
}

Matrix& Matrix::operator/=(const double d) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m[i][j] /= d;
    return *this;
}

Matrix operator - (const Matrix& mat) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = -mat[i][j];
    return temp;
}

Matrix operator * (const Matrix& mat, double d) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = mat[i][j] * d;
    return temp;
}

Matrix operator * (double d, const Matrix& mat) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = d * mat[i][j];
    return temp;
}

Matrix operator / (const Matrix& mat, double d) {
    int rows = mat.numRows();
    int cols = mat.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = mat[i][j] / d;
    return temp;
}

Matrix operator + (const Matrix& m1, const Matrix& m2) {
    int rows = m1.numRows();
    int cols = m1.numCols();
    if (rows != m2.numRows() || cols != m2.numCols())
        error("matrix dimension mismatch");
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = m1[i][j] + m2[i][j];
    return temp;
}

Matrix operator - (const Matrix& m1, const Matrix& m2) {
    int rows = m1.numRows();
    int cols = m1.numCols();
    if (rows != m2.numRows() || cols != m2.numCols())
        error("matrix dimension mismatch");
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[i][j] = m1[i][j] - m2[i][j];
    return temp;
}

Matrix operator * (const Matrix& m1, const Matrix& m2) {
    int n = m1.numCols();
    if (n != m2.numRows())
        error("matrix dimension mismatch");
    int rows = m1.numRows();
    int cols = m2.numCols();
    Matrix temp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            for (int k = 0; k < n; k++)
                temp[i][j] += m1[i][k] * m2[k][j];
    return temp;
}

Matrix Matrix::transpose() {
    Matrix temp(cols, rows);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            temp[j][i] = m[i][j];
    return temp;
}

std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            os << '\t' << mat.m[i][j];
        }
        os << '\n';
    }
    return os;
}

void solveGaussJordan(Matrix& A, Matrix& B) {

    int n = A.numRows();
    if (A.numCols() != n)
        error("matrix dimension mismatch");
    if (B.numRows() != n)
        error("matrix dimension mismatch");
    int m = B.numCols();

    int *indxc = new int [n];
    int *indxr = new int [n];
    int *ipiv = new int [n];

    for (int j = 0; j < n; j++)
        ipiv[j] = 0;

    for (int i = 0; i < n; i++) {
        double big = 0;
        int irow, icol;
        for (int j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        double Ajk = A[j][k];
                        if (std::abs(Ajk) >= big) {
                            big = std::abs(Ajk);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++ipiv[icol];

        if (irow != icol) {
            for (int j = 0; j < n; j++)
                std::swap(A[irow][j], A[icol][j]);
            for (int j = 0; j < m; j++)
                std::swap(B[irow][j], B[icol][j]);
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (A[icol][icol] == 0.0)
            error("solveGaussJordan singular matrix");
        double pivinv = 1 / A[icol][icol];
        A[icol][icol] = 1;
        for (int j = 0; j < n; j++)
            A[icol][j] *= pivinv;
        for (int j = 0; j < m; j++)
            B[icol][j] *= pivinv;

        for (int j = 0; j < n; j++) {
            if (j != icol) {
                double mult = A[j][icol];
                for (int k = 0; k < n; k++)
                    A[j][k] -= A[icol][k] * mult;
                for (int k = 0; k < m; k++)
                    B[j][k] -= B[icol][k] * mult;
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        if (indxr[i] != indxc[i]) {
            for (int j = 0; j < n; j++)
                std::swap(A[j][indxr[i]], A[j][indxc[i]]);
        }
    }

    delete [] indxc;
    delete [] indxr;
    delete [] ipiv;
}

void ludcmp(Matrix& A, int *indx, double& d) {
    const double TINY=1.0e-20;
    int i,imax,j,k;
    double big,dum,sum,temp;

    int n=A.numRows();
    Vector vv(n);
    d=1.0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if ((temp=std::abs(A[i][j])) > big) big=temp;
        if (big == 0.0) error("Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            sum=A[i][j];
            for (k=0;k<i;k++) sum -= A[i][k]*A[k][j];
            A[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            sum=A[i][j];
            for (k=0;k<j;k++) sum -= A[i][k]*A[k][j];
            A[i][j]=sum;
            if ((dum=vv[i]*std::abs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0;k<n;k++) {
                dum=A[imax][k];
                A[imax][k]=A[j][k];
                A[j][k]=dum;
            }
            d = -d;
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j][j] == 0.0) A[j][j]=TINY;
        if (j != n-1) {
            dum=1.0/(A[j][j]);
            for (i=j+1;i<n;i++) A[i][j] *= dum;
        }
    }
}

void lubksb(const Matrix& A, int *indx, Matrix& B) {
    int i,ii=0,ip,j;
    double sum;

    int n=A.numRows();
    for (int k = 0; k < B.numCols(); k++) {
        for (i=0;i<n;i++) {
            ip=indx[i];
            sum=B[ip][k];
            B[ip][k]=B[i][k];
            if (ii != 0)
                for (j=ii-1;j<i;j++) sum -= A[i][j]*B[j][k];
            else if (sum != 0.0)
                ii=i+1;
            B[i][k]=sum;
        }
        for (i=n-1;i>=0;i--) {
            sum=B[i][k];
            for (j=i+1;j<n;j++) sum -= A[i][j]*B[j][k];
            B[i][k]=sum/A[i][i];
        }
    }
}

void solveLUDecompose(Matrix& A, Matrix& B) {

    int n = A.numRows();
    if (A.numCols() != n)
        error("matrix dimension mismatch");
    if (B.numRows() != n)
        error("matrix dimension mismatch");
    int m = B.numCols();

    int *indx = new int [n];
    double d;
    ludcmp(A, indx, d);
    lubksb(A, indx, B);

    delete [] indx;
}

void reduceHouseholder(Matrix& A, Vector& d, Vector& e) {

    int n = d.dimension();

    for (int i = n - 1; i > 0; i--) {
        int l = i - 1;
        double h = 0;
        double scale = 0;
        if (l > 0) {
            for (int k = 0; k <= l; k++)
                scale += std::abs(A[i][k]);
            if (scale == 0.0)
                e[i] = A[i][l];
            else {
                for (int k = 0; k <= l; k++) {
                    A[i][k] /= scale;
                    h += A[i][k] * A[i][k];
                }
                double f = A[i][l];
                double g = (f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                A[i][l] = f - g;
                f = 0.0;
                for (int j = 0; j <= l; j++) {
                    A[j][i] = A[i][j] / h;
                    g = 0.0;
                    for (int k = 0; k <= j; k++)
                       g += A[j][k] * A[i][k];
                    for (int k = j + 1; k <= l; k++)
                       g += A[k][j] * A[i][k];
                    e[j] = g / h;
                    f += e[j] * A[i][j];
                }
                double hh = f / (h + h);
                for (int j = 0; j <= l; j++) {
                    f = A[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; k++)
                        A[j][k] -= f * e[k] + g * A[i][k];
                }
            }
        } else
            e[i] = A[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0]=0.0;
    for (int i = 0; i < n; i++) {
        if (d[i] != 0.0) {
            for (int j = 0; j < i; j++) {
                double g = 0;
                for (int k = 0; k < i; k++)
                    g += A[i][k] * A[k][j];
                for (int k = 0; k < i; k++)
                    A[k][j] -= g * A[k][i];
            }
        }
        d[i] = A[i][i];
        A[i][i] = 1.0;
        for (int j = 0; j < i; j++)
            A[j][i] = A[i][j] = 0.0;
    }
}

static double pythag(double a, double b) {
    double absa = std::abs(a);
    double absb = std::abs(b);
    if (absa > absb) {
        double ratio = absb / absa;
        return absa * std::sqrt(1 + ratio * ratio);
    } else {
        if (absb == 0.0)
            return 0.0;
        else {
            double ratio = absa / absb;
            return absb * std::sqrt(1 + ratio * ratio);
        }
    }
}

inline double sign(double a, double b) {
    if (b >= 0.0) {
        if (a >= 0.0)
            return a;
        else
            return -a;
    } else {
        if (a >= 0.0)
            return -a;
        else
            return a;
    }
}

void solveTQLI(Vector& d, Vector& e, Matrix& z) {

    int n = d.dimension();
    for (int i = 1; i < n; i++)
        e[i-1] = e[i];
    e[n-1] = 0.0;
    for (int l = 0; l < n; l++) {
        int iter = 0, m;
        do {
            for (m = l ; m < n-1; m++) {
                double dd = std::abs(d[m]) + std::abs(d[m+1]);
                if ((std::abs(e[m]) + dd) == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 30)
                     error("Too many iterations in solveTQLI");
                double g = (d[l+1] - d[l]) / (2.0 * e[l]);
                double r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + sign(r, g));
                double s = 1.0;
                double c = 1.0;
                double p = 0.0;
                int i;
                for (i = m-1; i >= l; i--) {
                    double f = s * e[i];
                    double b = c * e[i];
                    e[i+1] = r = pythag(f, g);
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i+1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i+1] = g + (p = s * r);
                    g = c * r - b;
                    for (int k = 0; k < n; k++) {
                        f = z[k][i+1];
                        z[k][i+1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l)
                    continue;
                 d[l] -= p;
                 e[l] = g;
                 e[m] = 0.0;
            }
        } while (m != l);
    }
}

void sortEigen(Vector& d, Matrix& z) {

    // sorts eigenvalues and eigenvector in descending order
    int n = d.dimension();
    if (z.numRows() != n || z.numCols() != n)
        error("Bad vector, matrix dimensions in sortEigen");
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        double p = d[k];
        for (int j = i; j < n; j++)
            if (d[j] >= p)
                p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
                p = z[j][i];
                z[j][i] = z[j][k];
                z[j][k] = p;
            }
        }
    }
}

void singularValueDecompose(Matrix& a, Vector& w, Matrix& v) {

    bool flag;
    int i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;

    int m=a.numRows();
    int n=a.numCols();
    Vector rv1(n);
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
        l=i+2;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += std::abs(a[k][i]);
            if (scale != 0.0) {
                for (k=i;k<m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -sign(std::sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i+1 <= m && i != n) {
            for (k=l-1;k<n;k++) scale += std::abs(a[i][k]);
            if (scale != 0.0) {
                for (k=l-1;k<n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l-1];
                g = -sign(std::sqrt(s),f);
                h=f*g-s;
                a[i][l-1]=f-g;
                for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
                for (j=l-1;j<m;j++) {
                    for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
                    for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l-1;k<n;k++) a[i][k] *= scale;
            }
        }
        anorm=std::max(anorm,(std::abs(w[i])+std::abs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g != 0.0) {
                for (j=l;j<n;j++)
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=std::min(m,n)-1;i>=0;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) a[i][j]=0.0;
        if (g != 0.0) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<m;j++) a[j][i] *= g;
        } else for (j=i;j<m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=true;
            for (l=k;l>=0;l--) {
                nm=l-1;
                if ((std::abs(rv1[l])+anorm) == anorm) {
                    flag=false;
                    break;
                }
                if ((std::abs(w[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l-1;i<k+1;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((std::abs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) error("singularValueDecomposition: no convergence in 30 iterations");
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z != 0.0) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
}

Vector solveEigenSymmetric(Matrix& A) {
    int n = A.numRows();
    Vector d(n), e(n);
    reduceHouseholder(A, d, e);
    solveTQLI(d, e, A);
    sortEigen(d, A);
    return d;
}

void lineSearch(Vector& xold, double fold, Vector& g, Vector& p, Vector& x,
                double& f, double stpmax, bool& check, double (*func)(Vector&))
{
    const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
    int i;
    double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    double rhs1,rhs2,slope,sum,temp,test,tmplam;

    int n = xold.dimension();
    check = false;
    sum=0.0;
    for (i=0;i<n;i++) sum += p[i]*p[i];
    sum=std::sqrt(sum);
    if (sum > stpmax)
        for (i=0;i<n;i++) p[i] *= stpmax/sum;
    slope = 0.0;
    for (i=0;i<n;i++)
        slope += g[i]*p[i];
    if (slope >= 0.0) error("Roundoff problem in lineSearch");
    test=0.0;
    for (i=0;i<n;i++) {
        temp=std::abs(p[i])/std::max(std::abs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
        f=func(x);
        if (alam < alamin) {
            for (i=0;i<n;i++) x[i]=xold[i];
            check = true;
            return;
        } else if (f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(f-fold-slope));
            else {
                rhs1 = f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc < 0.0) tmplam=0.5*alam;
                    else if (b <= 0.0) tmplam=(-b+std::sqrt(disc))/(3.0*a);
                    else tmplam=-slope/(b+std::sqrt(disc));
                }
                if (tmplam > 0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = f;
        alam=std::max(tmplam,0.1*alam);
    }
}

void minimizeBFGS(Vector& p, const double gtol, int& iter, double& fret,
            double (*func)(Vector&), void (*dfunc)(Vector&, Vector&)) {
    const int ITMAX = 200;
    const double EPS = std::numeric_limits<double>::epsilon();
    const double TOLX = 4*EPS, STPMX = 100.0;
    bool check;
    int i,its,j;
    double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

    int n = p.dimension();
    Vector dg(n),g(n),hdg(n),pnew(n),xi(n);
    Matrix hessin(n, n);

    fp=func(p);
    dfunc(p,g);
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) hessin[i][j]=0.0;
        hessin[i][i]=1.0;
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    stpmax=STPMX*std::max(std::sqrt(sum),double(n));
    for (its=0;its<ITMAX;its++) {
        iter=its;
        lineSearch(p,fp,g,xi,pnew,fret,stpmax,check,func);
        fp = fret;
        for (i=0;i<n;i++) {
            xi[i]=pnew[i]-p[i];
            p[i]=pnew[i];
        }
        test=0.0;
        for (i=0;i<n;i++) {
            temp=std::abs(xi[i])/std::max(std::abs(p[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX)
            return;
        for (i=0;i<n;i++) dg[i]=g[i];
        dfunc(p,g);
        test=0.0;
        den=std::max(fret,1.0);
        for (i=0;i<n;i++) {
            temp=std::abs(g[i])*std::max(std::abs(p[i]),1.0)/den;
            if (temp > test) test=temp;
        }
        if (test < gtol)
            return;
        for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
        for (i=0;i<n;i++) {
            hdg[i]=0.0;
            for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
        }
        fac=fae=sumdg=sumxi=0.0;
        for (i=0;i<n;i++) {
            fac += dg[i]*xi[i];
            fae += dg[i]*hdg[i];
            sumdg += dg[i]*dg[i];
            sumxi += xi[i]*xi[i];
        }
        if (fac > std::sqrt(EPS*sumdg*sumxi)) {
            fac=1.0/fac;
            fad=1.0/fae;
            for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
            for (i=0;i<n;i++) {
                for (j=i;j<n;j++) {
                    hessin[i][j] += fac*xi[i]*xi[j]
                        -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                    hessin[j][i]=hessin[i][j];
                }
            }
        }
        for (i=0;i<n;i++) {
            xi[i]=0.0;
            for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
        }
    }
    error("Too many iterations in minimizeBFGS");
}

} /* end namespace cpl */
