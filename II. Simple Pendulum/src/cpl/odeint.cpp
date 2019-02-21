#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "vector.hpp"
#include "odeint.hpp"

namespace cpl {

//  Fourth order Runge-Kutta
void RK4Step(Vector& x, double tau,
             Vector derivs(const Vector&))
{
    Vector k1 = tau * derivs(x);
    Vector k2 = tau * derivs(x + 0.5 * k1);
    Vector k3 = tau * derivs(x + 0.5 * k2);
    Vector k4 = tau * derivs(x + k3);
    x += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

//  Adaptive step size control using Runge-Kutta and step doubling
void adaptiveRK4Step(Vector& x, double& tau, double accuracy,
                     Vector derivs(const Vector&))
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    int n = x.dimension();
    Vector x_half(n), x_full(n), Delta(n);
    Vector scale = derivs(x);
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take two half steps
        double tau_half = tau / 2;
        x_half = x;
        RK4Step(x_half, tau_half, derivs);
        RK4Step(x_half, tau_half, derivs);
        // take full step
        x_full = x;
        RK4Step(x_full, tau, derivs);
        // estimate error
        Delta = x_half - x_full;
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = max(err_max, abs(Delta[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = max(tau_temp, 0.1 * tau);
        else
            tau = min(tau_temp, 0.1 * tau);
        if (abs(tau) == 0.0) {
            cerr << "adaptiveRK4Step: step size underflow\naborting ..."
                 << endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    x = x_half + Delta / 15.0;
}

//  Runge-Kutta-Cash-Karp including error estimate
static void rkck(Vector& x, double tau,
                 Vector derivs(const Vector&), Vector& x_err)
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

    Vector
        k1 = tau * derivs(x),
        k2 = tau * derivs(x + b21*k1),
        k3 = tau * derivs(x + b31*k1 + b32*k2),
        k4 = tau * derivs(x + b41*k1 + b42*k2 + b43*k3),
        k5 = tau * derivs(x + b51*k1 + b52*k2 + b53*k3 + b54*k4),
        k6 = tau * derivs(x + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5);
    x += c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6;
    x_err = dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;
}

//  Runge-Kutta-Cash-Karp step
void RKCKStep(Vector& x, double tau, Vector derivs(const Vector&))
{
    Vector x_err(x.dimension());
    rkck(x, tau, derivs, x_err);
}

//  Adaptive step size control using Runge-Kutta-Cash-Karp
void adaptiveRKCKStep(Vector& x, double& tau, double accuracy,
                      Vector derivs(const Vector&))
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    int n = x.dimension();
    Vector x_err(n), x_temp(n);
    Vector scale = derivs(x);
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take Cash-Karp step including error estimate
        x_temp = x;
        rkck(x_temp, tau, derivs, x_err);
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = max(err_max, abs(x_err[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = max(tau_temp, 0.1 * tau);
        else
            tau = min(tau_temp, 0.1 * tau);
        if (abs(tau) == 0.0) {
            cerr << "adaptiveRKCKStep: step size underflow\naborting ..."
                 << endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    x = x_temp;
}

} /* end namespace cpl */
