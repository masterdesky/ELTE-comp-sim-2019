#ifndef CPL_RANDOM_HPP
#define CPL_RANDOM_HPP

namespace cpl {

//  Random integer seed based on the current time of day
extern int timeSeed();

//  Standard C++ Library generator
extern double rand(int& seed, bool set=false);

//  Park-Miller generator with Schrage's algorithm to prevent overflow
extern double ranf(int& seed, bool set=false);

//  Quick and dirty generator according to Numerical Recipes
extern double qand(int& seed, bool set=false);

//  Gaussian distributed random number based on Box-Muller algorithm
extern double gasdev(int& seed);

}  /* end namespace cpl */

#endif /* CPL_RANDOM_HPP */
