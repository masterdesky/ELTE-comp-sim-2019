#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "optimize.hpp"

static    int sign_of_func;
static    double (*user_func)(double);
static    double func(double x) {
    return sign_of_func * (*user_func)(x);
}

inline void shft3(double &a, double &b, double &c, const double d) {
    a=b; b=c; c=d;
}

inline double SIGN(const double& a, const double& b) {
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

namespace cpl {

double Optimizer::findMaximum(const double a, const double b,
    double f(double), double tol, double& xmin) {
    sign_of_func = -1;
    user_func = f;
    double ax = a, bx = b, cx, fa, fb, fc;
    mnbrak(ax, bx, cx, fa, fb, fc, func);
    return -golden(ax, bx, cx, func, tol, xmin);
}

double Optimizer::findMinimum(const double a, const double b,
    double f(double), double tol, double& xmin) {
    double ax = a, bx = b, cx, fa, fb, fc;
    mnbrak(ax, bx, cx, fa, fb, fc, f);
    return golden(ax, bx, cx, f, tol, xmin);
}

void Optimizer::mnbrak(double &ax, double &bx, double &cx,
                       double &fa, double &fb, double &fc,
                       double func(const double)) {
    const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
    double ulim,u,r,q,fu;

    fa=func(ax);
    fb=func(bx);
    if (fb > fa) {
        swap(ax,bx);
        swap(fb,fa);
    }
    cx=bx+GOLD*(bx-ax);
    fc=func(cx);
    while (fb > fc) {
        r=(bx-ax)*(fb-fc);
        q=(bx-cx)*(fb-fa);
        u=bx-((bx-cx)*q-(bx-ax)*r)/
            (2.0*SIGN(max(abs(q-r),TINY),q-r));
        ulim=bx+GLIMIT*(cx-bx);
        if ((bx-u)*(u-cx) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                ax=bx;
                bx=u;
                fa=fb;
                fb=fu;
                return;
            } else if (fu > fb) {
                cx=u;
                fc=fu;
                return;
            }
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        } else if ((cx-u)*(u-ulim) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                shft3(bx,cx,u,cx+GOLD*(cx-bx));
                shft3(fb,fc,fu,func(u));
            }
        } else if ((u-ulim)*(ulim-cx) >= 0.0) {
            u=ulim;
            fu=func(u);
        } else {
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        }
        shft3(ax,bx,cx,u);
        shft3(fa,fb,fc,fu);
    }

}

inline void shft2(double &a, double &b, const double c) {
    a=b;  b=c;
}

double Optimizer::golden(const double ax, const double bx, const double cx,
                         double f(const double), const double tol,
                         double &xmin) {
    const double R=0.61803399,C=1.0-R;
    double f1,f2,x0,x1,x2,x3;

    x0=ax;
    x3=cx;
    if (abs(cx-bx) > abs(bx-ax)) {
        x1=bx;
        x2=bx+C*(cx-bx);
    } else {
        x2=bx;
        x1=bx-C*(bx-ax);
    }
    f1=f(x1);
    f2=f(x2);
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
        if (f2 < f1) {
            shft3(x0,x1,x2,R*x2+C*x3);
            shft2(f1,f2,f(x2));
        } else {
            shft3(x3,x2,x1,R*x1+C*x0);
            shft2(f2,f1,f(x1));
        }
    }
    if (f1 < f2) {
        xmin=x1;
        return f1;
    } else {
        xmin=x2;
        return f2;
    }
}

} /* end namespace cpl */
