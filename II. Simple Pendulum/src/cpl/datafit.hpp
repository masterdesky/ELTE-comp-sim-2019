#ifndef CPL_DATAFIT_HPP
#define CPL_DATAFIT_HPP

namespace cpl {

class Vector;

extern void fit(          // Least Squares fit to y(x) = a + b x
    const Vector& x,      // Input: Vector of x values
    const Vector& y,      // Input: Vector of y values
    const Vector& sig,    // Input: Vector of individual standard deviations
    const bool mwt,       // Input: if false sig's assumed to be unavailable
    double& a,            // Output: intercept a
    double& b,            // Output: slope b
    double& siga,         // Output: uncertainty in a
    double& sigb,         // Output: uncertainty in b
    double& chi2,         // Output: chi-square
    double& q             // Output: goodness of fit
);

}  /* end namespace cpl */

#endif /* CPL_DATAFIT_HPP */
