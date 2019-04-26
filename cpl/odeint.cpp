#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vector.hpp"
#include "odeint.hpp"

namespace cpl {

// #################### SIMPLE ####################

//  Simple Euler
void EulerStep(cpl::Vector& x, double tau,
               cpl::Vector derivates(const cpl::Vector&))
{   
    // Time
    // Deflection
    // Velocity
    x += tau * derivates(x);
}

//  Semi-Implicit, Euler-Cromer
void EulerCromerStep(cpl::Vector& x, double tau,
                     cpl::Vector derivates(const cpl::Vector&))
{
    auto x_temp = x;

    // Time
    x_temp[0] += tau * derivates(x)[0];

    // Velocity
    x_temp[2] += tau * derivates(x)[2];

    // Deflection
    x_temp[1] += tau * x_temp[2];

    x = x_temp;
}

//  Fourth order Runge-Kutta
void RK4Step(cpl::Vector& x, double tau,
             cpl::Vector derivates(const cpl::Vector&))
{
    cpl::Vector 
        k1 = tau * derivates(x),
        k2 = tau * derivates(x + 0.5 * k1),
        k3 = tau * derivates(x + 0.5 * k2),
        k4 = tau * derivates(x + k3);

    x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

//  Adaptive step size control using Runge-Kutta and step doubling
void adaptiveRK4Step(cpl::Vector& x, double& tau, double accuracy,
                     cpl::Vector derivates(const cpl::Vector&))
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    int n = x.dimension();
    cpl::Vector x_half(n), x_full(n), Delta(n);
    cpl::Vector scale = derivates(x);
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take two half steps
        double tau_half = tau / 2;
        x_half = x;
        RK4Step(x_half, tau_half, derivates);
        RK4Step(x_half, tau_half, derivates);
        // take full step
        x_full = x;
        RK4Step(x_full, tau, derivates);
        // estimate error
        Delta = x_half - x_full;
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = std::max(err_max, abs(Delta[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = std::max(tau_temp, 0.1 * tau);
        else
            tau = std::min(tau_temp, 0.1 * tau);
        if (abs(tau) == 0.0) {
            std::cerr << "adaptiveRK4Step: step size underflow\naborting ..."
                 << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    x = x_half + Delta / 15.0;
}

//  Runge-Kutta-Cash-Karp including error estimate
static void rkck(cpl::Vector& x, double tau,
                 cpl::Vector derivates(const cpl::Vector&), cpl::Vector& x_err)
{
    const double b21 = 1.0/5.0, b31 = 3.0/40.0, b41 = 3.0/10.0,
        b51 = -11.0/54.0, b61 = 1631.0/55296.0, b32 = 9.0/40.0,
        b42 = -9.0/10.0, b52 = 5.0/2.0, b62 = 175.0/512.0, b43 = 6.0/5.0,
        b53 = -70.0/27.0, b63 = 575.0/13824.0, b54 = 35.0/27.0,
        b64 = 44275.0/110592.0, b65 = 253.0/4096.0, c1 = 37.0/378.0,
        c2 = 0.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c5 = 0.0,
        c6 = 512.0/1771.0, dc1 = c1 - 2825.0/27648.0, dc2 = c2 - 0.0,
        dc3 = c3 - 18575.0/48384.0, dc4 = c4 - 13525.0/55296.0,
        dc5 = c5 - 277.0/14336.0, dc6 = c6 - 1.0/4.0;

    cpl::Vector
        k1 = tau * derivates(x),
        k2 = tau * derivates(x + b21*k1),
        k3 = tau * derivates(x + b31*k1 + b32*k2),
        k4 = tau * derivates(x + b41*k1 + b42*k2 + b43*k3),
        k5 = tau * derivates(x + b51*k1 + b52*k2 + b53*k3 + b54*k4),
        k6 = tau * derivates(x + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5);
    x += c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6;
    x_err = dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;
}

//  Runge-Kutta-Cash-Karp step
void RKCKStep(cpl::Vector& x, double tau,
              cpl::Vector derivates(const cpl::Vector&))
{
    cpl::Vector x_err(x.dimension());
    rkck(x, tau, derivates, x_err);
}

//  Adaptive step size control using Runge-Kutta-Cash-Karp
void adaptiveRKCKStep(cpl::Vector& x, double& tau, double accuracy,
                      cpl::Vector derivates(const cpl::Vector&))
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    int n = x.dimension();
    cpl::Vector x_err(n), x_temp(n);
    cpl::Vector scale = derivates(x);
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take Cash-Karp step including error estimate
        x_temp = x;
        rkck(x_temp, tau, derivates, x_err);
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = std::max(err_max, abs(x_err[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = std::max(tau_temp, 0.1 * tau);
        else
            tau = std::min(tau_temp, 0.1 * tau);
        if (abs(tau) == 0.0) {
            std::cerr << "adaptiveRKCKStep: step size underflow\naborting ..."
                 << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    x = x_temp;
}

// #################### DOUBLE ####################

//  Simple Euler
void EulerStep_double(cpl::Vector& x, cpl::Vector& y, double tau,
                      std::vector<cpl::Vector> derivates(const cpl::Vector&, const cpl::Vector&))
{   
    cpl::Vector x_temp = x;
    // Time_1, Deflection_1, Velocity_1
    x += tau * (derivates(x, y)[0]);

    // Time_2, Deflection_2, Velocity_2
    y += tau * (derivates(x_temp, y)[1]);
}


//  Semi-Implicit, Euler-Cromer
void EulerCromerStep_double(cpl::Vector& x, cpl::Vector& y, double tau,
                            std::vector<cpl::Vector> derivates(const cpl::Vector&, const cpl::Vector&))
{
    cpl::Vector x_temp = x;
    // Time
    x[0] += tau * (derivates(x, y)[0])[0];
    // Velocity
    x[2] += tau * (derivates(x, y)[0])[2];
    // Deflection
    x[1] += tau * x[2];

    // Time
    y[0] += tau * (derivates(x_temp, y)[1])[0];
    // Velocity
    y[2] += tau * (derivates(x_temp, y)[1])[2];
    // Deflection
    y[1] += tau * y[2];
}


//  Fourth order Runge-Kutta
void RK4Step_double(cpl::Vector& x, cpl::Vector& y, double tau,
                    std::vector<cpl::Vector> derivates(const cpl::Vector&, const cpl::Vector&))
{
    std::vector<cpl::Vector>
        k1 = {tau * derivates(x, y)[0], tau * derivates(x, y)[1]},
        k2 = {tau * derivates(x + 0.5 * k1[0], y + 0.5 * k1[1])[0], tau * derivates(x + 0.5 * k1[0], y + 0.5 * k1[1])[1]},
        k3 = {tau * derivates(x + 0.5 * k2[0], y + 0.5 * k2[1])[0], tau * derivates(x + 0.5 * k2[0], y + 0.5 * k2[1])[1]},
        k4 = {tau * derivates(x + k3[0], y + k3[1])[0]            , tau * derivates(x + k3[0], y + k3[1])[1]};
    x += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0;
    y += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0;
}


//  Runge-Kutta-Cash-Karp including error estimate
static void rkck_double(cpl::Vector& x, cpl::Vector& y, double tau,
                        std::vector<cpl::Vector> derivates(const cpl::Vector&, const cpl::Vector&),
                        cpl::Vector& x_err, cpl::Vector& y_err)
{
    const double b21 = 1.0/5.0, b31 = 3.0/40.0, b41 = 3.0/10.0,
        b51 = -11.0/54.0, b61 = 1631.0/55296.0, b32 = 9.0/40.0,
        b42 = -9.0/10.0, b52 = 5.0/2.0, b62 = 175.0/512.0, b43 = 6.0/5.0,
        b53 = -70.0/27.0, b63 = 575.0/13824.0, b54 = 35.0/27.0,
        b64 = 44275.0/110592.0, b65 = 253.0/4096.0, c1 = 37.0/378.0,
        c2 = 0.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c5 = 0.0,
        c6 = 512.0/1771.0, dc1 = c1 - 2825.0/27648.0, dc2 = c2 - 0.0,
        dc3 = c3 - 18575.0/48384.0, dc4 = c4 - 13525.0/55296.0,
        dc5 = c5 - 277.0/14336.0, dc6 = c6 - 1.0/4.0;

    std::vector<cpl::Vector>
        k1 = {tau * derivates(x, y)[0],
              tau * derivates(x, y)[1]},

        k2 = {tau * derivates(x + b21*k1[0], y + b21*k1[1])[0],
              tau * derivates(x + b21*k1[0], y + b21*k1[1])[1]},

        k3 = {tau * derivates(x + b31*k1[0] + b32*k2[0], y + b31*k1[1] + b32*k2[1])[0],
              tau * derivates(x + b31*k1[0] + b32*k2[0], y + b31*k1[1] + b32*k2[1])[1]},

        k4 = {tau * derivates(x + b41*k1[0] + b42*k2[0] + b43*k3[0], y + b41*k1[1] + b42*k2[1] + b43*k3[1])[0],
              tau * derivates(x + b41*k1[0] + b42*k2[0] + b43*k3[0], y + b41*k1[1] + b42*k2[1] + b43*k3[1])[1]},

        k5 = {tau * derivates(x + b51*k1[0] + b52*k2[0] + b53*k3[0] + b54*k4[0], y + b51*k1[1] + b52*k2[1] + b53*k3[1] + b54*k4[1])[0],
              tau * derivates(x + b51*k1[0] + b52*k2[0] + b53*k3[0] + b54*k4[0], y + b51*k1[1] + b52*k2[1] + b53*k3[1] + b54*k4[1])[1]},

        k6 = {tau * derivates(x + b61*k1[0] + b62*k2[0] + b63*k3[0] + b64*k4[0] + b65*k5[0], y + b61*k1[1] + b62*k2[1] + b63*k3[1] + b64*k4[1] + b65*k5[1])[0],
              tau * derivates(x + b61*k1[0] + b62*k2[0] + b63*k3[0] + b64*k4[0] + b65*k5[0], y + b61*k1[1] + b62*k2[1] + b63*k3[1] + b64*k4[1] + b65*k5[1])[1]};

    x += c1*k1[0] + c2*k2[0] + c3*k3[0] + c4*k4[0] + c5*k5[0] + c6*k6[0];
    y += c1*k1[1] + c2*k2[1] + c3*k3[1] + c4*k4[1] + c5*k5[1] + c6*k6[1];
    x_err = dc1*k1[0] + dc2*k2[0] + dc3*k3[0] + dc4*k4[0] + dc5*k5[0] + dc6*k6[0];
    y_err = dc1*k1[1] + dc2*k2[1] + dc3*k3[1] + dc4*k4[1] + dc5*k5[1] + dc6*k6[1];
}


//  Runge-Kutta-Cash-Karp step
void RKCKStep_double(cpl::Vector& x, cpl::Vector& y, double tau,
                     std::vector<cpl::Vector> derivates(const cpl::Vector&, const cpl::Vector&))
{
    cpl::Vector x_err(x.dimension());
    cpl::Vector y_err(y.dimension());
    rkck_double(x, y, tau, derivates, x_err, y_err);
}

} /* end namespace cpl */
