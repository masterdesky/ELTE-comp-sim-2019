#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

#include "interpolate.hpp"

namespace cpl {

// polynomial interpolation routine adapted from Numerical Recipes
// polynomial degree inferred from dimension of input vector

void Interpolator::polint(        // uses Neville's algorithm
    const vector<double>& xa,     // vector of abscissa values
    const vector<double>& ya,     // vector of ordinate values
    const double x,               // point x
    double& y,                    // interpolated y(x)
    double& dy)                   // estimate of error in interpolation
{
    int i,m,ns=0;
    double den,dif,dift,ho,hp,w;

    int n=xa.size();
    vector<double> c(n),d(n);
    dif=abs(x-xa[0]);
    for (i=0;i<n;i++) {
        if ((dift=abs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=0;i<n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ((den=ho-hp) == 0.0) { // nrerror("Error in routine polint");
                cerr << "Error in Interpolator::polint\n"
                     << "Neville's algorithm breaks down when two "
                     << "input x values are equal within roundoff\n"
                     << "x, y values" << endl;
                for (int j = 0; j < n; j++)
                    cerr << xa[j] << '\t' << ya[j] << '\n';
                exit(EXIT_FAILURE);
            }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
    }
}

}  /* end namespace cpl */
