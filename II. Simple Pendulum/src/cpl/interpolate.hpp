#ifndef CPL_INTERPOLATE_HPP
#define CPL_INTERPOLATE_HPP

namespace cpl {

class Interpolator {
  public:

    static void polint(             // adapted from Numerical Recipes
        const vector<double>& xa,   // vector of x values - input
        const vector<double>& ya,   // vector of y values - input
        const double x,             // x at which to find y - input
        double& y,                  // value of y - output
        double& dy);                // estimate or error in y
};

}  /* end namespace cpl */

#endif /* CPL_INTERPOLATE_HPP */
