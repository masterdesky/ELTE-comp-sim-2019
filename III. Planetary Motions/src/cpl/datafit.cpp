#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;

#include "vector.hpp"
#include "datafit.hpp"

namespace cpl {

inline double SQR(double a) {
    return a * a;
}

static void error(const char *msg) {
    cerr << msg << endl;
    exit(EXIT_FAILURE);
}

double gammln(const double xx)
{
    int j;
    double x,y,tmp,ser;
    static const double cof[6]=
        {76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
         -0.5395239384953e-5};

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<6;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void gser(double& gamser, const double a, const double x, double& gln)
{
    const int ITMAX=100;
    const double EPS=numeric_limits<double>::epsilon();
    int n;
    double sum,del,ap;

    gln=gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) error("x less than 0 in routine gser");
        gamser=0.0;
        return;
    } else {
        ap=a;
        del=sum=1.0/a;
        for (n=0;n<ITMAX;n++) {
            ++ap;
            del *= x/ap;
            sum += del;
            if (abs(del) < abs(sum)*EPS) {
                gamser=sum*exp(-x+a*log(x)-gln);
                return;
            }
        }
        error("a too large, ITMAX too small in routine gser");
        return;
    }
}

void gcf(double& gammcf, const double a, const double x, double& gln)
{
    const int ITMAX=100;
    const double EPS=numeric_limits<double>::epsilon();
    const double FPMIN=numeric_limits<double>::min()/EPS;
    int i;
    double an,b,c,d,del,h;

    gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
        an = -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (abs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (abs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (abs(del-1.0) <= EPS) break;
    }
    if (i > ITMAX) error("a too large, ITMAX too small in gcf");
    gammcf=exp(-x+a*log(x)-gln)*h;
}

double gammq(const double a, const double x)
{
    double gamser,gammcf,gln;

    if (x < 0.0 || a <= 0.0)
        error("Invalid arguments in routine gammq");
    if (x < a+1.0) {
        gser(gamser,a,x,gln);
        return 1.0-gamser;
    } else {
        gcf(gammcf,a,x,gln);
        return gammcf;
    }
}

void fit(
    const Vector& x,
    const Vector& y,
    const Vector& sig,
    const bool mwt,
    double& a,
    double& b,
    double& siga,
    double& sigb,
    double& chi2,
    double& q)
{
    int i;
    double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

    int ndata=x.size();
    b=0.0;
    if (mwt) {
        ss=0.0;
        for (i=0;i<ndata;i++) {
            wt=1.0/SQR(sig[i]);
            ss += wt;
            sx += x[i]*wt;
            sy += y[i]*wt;
        }
    } else {
        for (i=0;i<ndata;i++) {
            sx += x[i];
            sy += y[i];
        }
        ss=ndata;
    }
    sxoss=sx/ss;
    if (mwt) {
        for (i=0;i<ndata;i++) {
            t=(x[i]-sxoss)/sig[i];
            st2 += t*t;
            b += t*y[i]/sig[i];
        }
    } else {
        for (i=0;i<ndata;i++) {
            t=x[i]-sxoss;
            st2 += t*t;
            b += t*y[i];
        }
    }
    b /= st2;
    a=(sy-sx*b)/ss;
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    sigb=sqrt(1.0/st2);
    chi2=0.0;
    q=1.0;
    if (!mwt) {
        for (i=0;i<ndata;i++)
            chi2 += SQR(y[i]-a-b*x[i]);
        sigdat=sqrt(chi2/(ndata-2));
        siga *= sigdat;
        sigb *= sigdat;
    } else {
        for (i=0;i<ndata;i++)
            chi2 += SQR((y[i]-a-b*x[i])/sig[i]);
        if (ndata>2) q=gammq(0.5*(ndata-2),0.5*chi2);
    }
}

}  /* end namespace cpl */
