#ifndef CPL_OPTIMIZE_HPP
#define CPL_OPTIMIZE_HPP

namespace cpl {

class Optimizer {

  public:

    double findMaximum(         // returns f(x) at local maximum
        const double xa,        // guess for one point near desried maximum
        const double xb,        // guess for second point near maximum
        double f(double),       // name of function to be maximized
        double accuracy,        // desired accuracy in x
        double& xMax);          // value of x at local maximum

    double findMinimum(         // returns f(x) at local minimum
        const double xa,        // guess for one point near desried minimum
        const double xb,        // guess for second point near minimum
        double f(double),       // name of function to be minimized
        double accuracy,        // desired accuracy in x
        double& xMin);          // value of x at local minimum

    double golden(              // numerical recipes routine
        const double,
        const double,
        const double,
        double f(double),
        const double,
        double&
    );

    void mnbrak(                // numerical recipes bracketing routine
        double&,
        double&,
        double&,
        double&,
        double&,
        double&,
        double f(double)
    );

};

} /* end namespace cpl */

#endif /* CPL_OPTIMIZE_HPP */
