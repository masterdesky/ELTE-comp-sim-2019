#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

double L = 1.0;              // length of pendulum
double q = 0.5;              // damping coefficient
double Omega_D = 2.0/3.0;    // frequency of driving force
double F_D = 0.9;            // amplitude of driving force
bool nonlinear;              // linear if false

cpl::Vector f(const cpl::Vector& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    cpl::Vector f(3);             // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (nonlinear)
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

int main(int argc, char* argv[]) {

    double theta, omega, t_max;

    std::cout << " Nonlinear damped driven pendulum\n"
              << " --------------------------------\n";
    std::string response = argv[1];
    nonlinear = (response[0] == 'n');

    L = atof(argv[2]);                  // Length of pendulum L
    q = atof(argv[3]);                  // Damping coefficient q
    Omega_D = atof(argv[4]);            // Driving frequencey Omega_D
    F_D = atof(argv[5]);                // Driving amplitude F_D
    theta = atof(argv[6]);              // Theta(0)
    omega = atof(argv[7]);              // Omega(ω)
    t_max = atof(argv[8]);              // Integration time t_max

    double dt = 0.05;
    double accuracy = 1e-6;
    std::ofstream dataFile("..\\out\\pendulum.dat");

    double t = 0;
    cpl::Vector x(3);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;

    while (t < t_max) {
        cpl::adaptiveRK4Step(x, dt, accuracy, f);
        t = x[0], theta = x[1], omega = x[2];
        if (nonlinear) {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        // Returns in radian (0) and radian\s (ω)
        dataFile << t << '\t' << theta << '\t' << omega << '\n';
    }

    std::cout << " Output data to file pendulum.dat" << std::endl;
    dataFile.close();
}

