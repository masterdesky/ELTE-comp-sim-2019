#ifndef CPL_ODEINT_HPP
#define CPL_ODEINT_HPP

#include <iostream>
#include <vector>

namespace cpl {

class Vector;

// ODE integration routines adapted from Numerical Recipes

// #################### SINGLE BODY ####################

//  take a single Euler step
extern void EulerStep(
    Vector& x,                        //  extended solution vector
    double tau,                       //  step size
    Vector derivs(const Vector& x)    //  extended derivative vector
);

//  take a single Euler step
extern void EulerCromerStep(
    Vector& x,                        //  extended solution vector
    double tau,                       //  step size
    Vector derivs(const Vector& x)    //  extended derivative vector
);

//  take a single 4th order Runge-Kutta step
extern void RK4Step(
    Vector& x,                        //  extended solution vector
    double tau,                       //  step size
    Vector derivs(const Vector& x)    //  extended derivative vector
);

//  take one adaptive RK4 step using step doubling
extern void adaptiveRK4Step(
    Vector& x,                        //  extended solution vector
    double& tau,                      //  recommended step size
    double accuracy,                  //  desired solution accuracy
    Vector derivs(const Vector&)      //  extended derivative vector
);

//  take one Runge-Kutta-Cash-Karp step
extern void RKCKStep(
    Vector& x,                        //  extended solution vector
    double tau,                       //  step size
    Vector derivs(const Vector& x)    //  derivative vector
);

//  take one adaptive RKCK step using embedded error estimate
extern void adaptiveRKCKStep(
    Vector& x,                        //  extended solution vector
    double& tau,                      //  recommended step size
    double accuracy,                  //  desired solution accuracy
    Vector derivs(const Vector&)      //  extended derivative vector
);

// #################### DOUBLE BODY ####################

//  take a single Euler step
extern void EulerStep_double(
    Vector& x,                        //  extended solution vector_1
    Vector& y,                        //  extended solution vector_2
    double tau,                       //  step size
    std::vector<cpl::Vector> derivs(const Vector& x,
                  const Vector& y)    //  extended derivative vector
);

//  take a single Euler step
extern void EulerCromerStep_double(
    Vector& x,                        //  extended solution vector_1
    Vector& y,                        //  extended solution vector_2
    double tau,                       //  step size
    std::vector<cpl::Vector> derivs(const Vector& x,
                  const Vector& y)    //  extended derivative vector
);

//  take a single 4th order Runge-Kutta step
extern void RK4Step_double(
    Vector& x,                        //  extended solution vector_1
    Vector& y,                        //  extended solution vector_2
    double tau,                       //  step size
    std::vector<cpl::Vector> derivs(const Vector& x,
                  const Vector& y)    //  extended derivative vector
);

//  take one Runge-Kutta-Cash-Karp step
extern void RKCKStep_double(
    Vector& x,                        //  extended solution vector_1
    Vector& y,                        //  extended solution vector_2
    double tau,                       //  step size
    std::vector<cpl::Vector> derivs(const Vector& x,
                  const Vector& y)    //  extended derivative vector
);

} /* end namespace cpl */

#endif /* CPL_ODEINT_HPP */
