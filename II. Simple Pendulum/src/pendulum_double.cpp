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

static double m_1;              // Mass of pendulum m_1
static double m_2;              // Mass of pendulum m_2
static double L_2;              // Length of pendulum L_2
static double L_1;              // Length of pendulum L_1
static double q_1;              // Damping coefficient q_1
static double q_2;              // Damping coefficient q_2
static double Omega_D_1;        // Driving frequency Omega_D_1
static double Omega_D_2;        // Driving frequency Omega_D_2
static double F_D_1;            // Driving amplitude F_D_1
static double F_D_2;            // Driving amplitude F_D_2
static double theta_1;          // Angle of deflection of pendulum (0_1)
static double theta_2;          // Angle of deflection of pendulum (0_2)
static double omega_1;          // Velocity of pendulum (ω_1)
static double omega_2;          // Velocity of pendulum (ω_2)

static double t;                // Current time t
static double t_max;            // Integration time t_max
static double dt;               // Stepsize
static double accuracy;         // Accuracy of simulation
static double error;            // Error of simulation


double energy_1() {
    return 0.5 * m_1 * (L_1 * L_1 * omega_1 * omega_1);
}

double energy_2() {
    return 0.5 * m_2 * (L_2 * L_2 * omega_2 * omega_2);
}

cpl::Vector f(const cpl::Vector& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    cpl::Vector f(3);                  // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (linearity)
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

int main(int argc, char* argv[]) {

    
    std::cout << " Nonlinear damped driven pendulum\n"
              << " --------------------------------\n";
    std::string response_linearity = argv[1];
    linearity = (response_linearity[0] == 'n');

    std::string mode = argv[2];

    L_2 = atof(argv[3]);                // Length of pendulum L_2
    L_1 = atof(argv[4]);                // Length of pendulum L_1
    q_1 = atof(argv[5]);                // Damping coefficient q_1
    q_2 = atof(argv[6]);                // Damping coefficient q_2
    Omega_D_1 = atof(argv[7]);          // Driving frequency Omega_D_1
    Omega_D_2 = atof(argv[8]);          // Driving frequency Omega_D_2
    F_D_1 = atof(argv[9]);              // Driving amplitude F_D_1
    F_D_2 = atof(argv[10]);             // Driving amplitude F_D_2
    theta_1 = atof(argv[11]);           // Angle of deflection of pendulum (0_1)
    theta_2 = atof(argv[12]);           // Angle of deflection of pendulum (0_2)
    omega_1 = atof(argv[13]);           // Velocity of pendulum (ω_1)
    omega_2 = atof(argv[14]);           // Velocity of pendulum (ω_2)

    t_max = atof(argv[15]);             // Integration time t_max
    dt = atof(argv[16]);                // Stepsize
    accuracy = atof(argv[17]);          // Accuracy of simulation
    
    std::ofstream dataFile("..\\out\\pendulum_double.dat");

    // Define containers for parameters of motion, and starting positions
    t = 0;
    cpl::Vector x(3);
    cpl::Vector y(3);

    x[0] = 0;
    x[1] = theta_1;
    x[2] = omega_1;

    y[0] = 0;
    y[1] = theta_2;
    y[2] = omega_2;

    // Time measurement init and starts
    auto start = std::chrono::steady_clock::now();
    std::chrono::microseconds duration;

    // STARTING ITERATION: t = 0
    // Returns in radian (0) and radian\s (ω)
    dataFile << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2 << '\t' << dt << '\t' << energy_1() << '\t' << energy_2() << '\t' << 0 << '\n';

    while (t < t_max) {
        if(mode == "runge") {
            cpl::RK4Step(x, dt, f);
        }
        else if(mode == "adapt_runge") {
            cpl::adaptiveRK4Step(x, dt, accuracy, f);
        }

        else if(mode == "rkck") {
            cpl::RKCKStep(x, dt, f);
        }
        else if(mode == "adapt_rkck") {
            cpl::adaptiveRKCKStep(x, dt, accuracy, f);
        }

        else if(mode == "euler") {
            cpl::EulerStep(x, dt, f);
        }

        else if(mode == "eulercromer") {
            cpl::EulerCromerStep(x, dt, f);
        }

        t = x[0];
        theta_1 = x[1], omega_1 = x[2];
        theta_2 = y[1], omega_2 = y[2];
        error = dt;

        if (linearity) {
            while (theta_1 >= pi) theta_1 -= 2 * pi;
            while (theta_1 < -pi) theta_1 += 2 * pi;
            while (theta_2 >= pi) theta_2 -= 2 * pi;
            while (theta_2 < -pi) theta_2 += 2 * pi;
        }
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
        // Returns in radian (0) and radian\s (ω)
        dataFile << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2 << '\t' << error << '\t' << energy_1() << '\t' << energy_2() << '\t' << duration.count() << '\n';
    }

    std::cout << " Output data to file pendulum_double.dat" << std::endl;
    dataFile.close();
}

