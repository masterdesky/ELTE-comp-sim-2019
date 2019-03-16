#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...

const double pi = 4 * atan(1.0);

const double g = 9.8;           // acceleration of gravity

bool linearity;                 // linear if false

static double m;                // Mass of pendulum m
static double L;                // Length of pendulum L
static double q;                // Damping coefficient q
static double Omega_D;          // Driving frequency Omega_D
static double F_D;              // Driving amplitude F_D
static double theta;            // Angle of deflection of pendulum (0)
static double omega;            // Velocity of pendulum (ω)

static double t;                // Current time t
static double t_max;            // Integration time t_max
static double dt;               // Stepsize
static double accuracy;         // Accuracy of simulation
static double error;            // Error of simulation


double energy() {
    return 0.5 * m * (L * L * omega * omega);
}

cpl::Vector derivates(const cpl::Vector& x) {  // extended derivative vector
    t = x[0];
    theta = x[1];
    omega = x[2];
    cpl::Vector f(3);             // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (linearity)
        f[2] = - m * (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - m * (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

int main(int argc, char* argv[]) {

    
    std::cout << " Nonlinear damped driven pendulum\n"
              << " --------------------------------\n";
    std::string response = argv[1];
    linearity = (response[0] == 'n');

    std::string mode = argv[2];

    m = atof(argv[3]);                  // Mass of pendulum m
    L = atof(argv[4]);                  // Length of pendulum L
    q = atof(argv[5]);                  // Damping coefficient q
    Omega_D = atof(argv[6]);            // Driving frequency Omega_D
    F_D = atof(argv[7]);                // Driving amplitude F_D
    theta = atof(argv[8]);              // Angle of deflection of pendulum (0)
    omega = atof(argv[9]);              // Velocity of pendulum (ω)
    t_max = atof(argv[10]);             // Integration time t_max
    dt = atof(argv[11]);                // Stepsize
    accuracy = atof(argv[12]);          // Accuracy of simulation
    
    std::ofstream dataFile("..\\out\\pendulum.dat");

    // Define containers for parameters of motion, and starting position
    t = 0;
    cpl::Vector x(3);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;

    // Time measurement init and starts
    auto start = std::chrono::steady_clock::now();
    std::chrono::microseconds duration;

    // STARTING ITERATION: t = 0
    // Returns in radian (0) and radian\s (ω)
    dataFile << t << '\t' << theta << '\t' << omega << '\t' << dt << '\t' << energy() << '\t' << 0 << '\n';

    while (t < t_max) {
        if(mode == "runge") {
            cpl::RK4Step(x, dt, derivates);
        }
        else if(mode == "adapt_runge") {
            cpl::adaptiveRK4Step(x, dt, accuracy, derivates);
        }

        else if(mode == "rkck") {
            cpl::RKCKStep(x, dt, derivates);
        }
        else if(mode == "adapt_rkck") {
            cpl::adaptiveRKCKStep(x, dt, accuracy, derivates);
        }

        else if(mode == "euler") {
            cpl::EulerStep(x, dt, derivates);
        }

        else if(mode == "eulercromer") {
            cpl::EulerCromerStep(x, dt, derivates);
        }

        t = x[0], theta = x[1], omega = x[2];
        error = dt;

        if (linearity) {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }

        // Runtime test
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        
        // Returns in radian (0) and radian\s (ω)
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << error << '\t' << energy() << '\t' << duration.count() << '\n';
    }

    std::cout << " Output data to file pendulum.dat" << std::endl;
    dataFile.close();
}

