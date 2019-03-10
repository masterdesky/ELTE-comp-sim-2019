#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...

const double pi = 4 * atan(1.0);

const double g = 9.8;           // acceleration of gravity

static double m_1;              // Mass of pendulum m_1
static double m_2;              // Mass of pendulum m_2
static double L_1;              // Length of pendulum L_1
static double L_2;              // Length of pendulum L_2
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
static double t_iteration;      // Integration time t_iteration
static double dt;               // Stepsize
static double accuracy;         // Accuracy of simulation
static double error;            // Error of simulation


double energy_1() {
    return 0.5 * m_1 * (L_1 * L_1 * omega_1 * omega_1);
}

double energy_2() {
    return 0.5 * m_2 * (L_2 * L_2 * omega_2 * omega_2);
}

std::vector<cpl::Vector> derivates(const cpl::Vector& x, const cpl::Vector& y) {             // extended derivative vectors
    t = x[0];
    theta_1 = x[1]; theta_2 = y[1];
    omega_1 = x[2]; omega_2 = y[2];
    cpl::Vector f_1(3); cpl::Vector f_2(3);                             // Vectors with 3 components
    
    f_1[0] = 1, f_2[0] = 1;
    f_1[1] = omega_1, f_2[1] = omega_2;

    // Expression for easier eq. of motions
    double d_theta = theta_1 - theta_2;

    double alfa = q_1 * omega_1 - F_D_1 * sin(Omega_D_1 * t);           // Damping + Driving for first pendulum
    double beta = q_2 * omega_2 - F_D_2 * sin(Omega_D_2 * t);           // Damping + Driving for second pendulum
    
    double gamma_1 = 2 * alfa - 2 * beta * cos(d_theta);
    double gamma_2 = 2 * alfa * cos(d_theta) - (2 * (m_1 + m_2)/m_2) * beta;

    // Eq. of motion for first pendulum
    f_1[2] = ((m_2 * L_1 * (omega_1 * omega_1) * sin(2 * d_theta) +
               2 * m_2 * L_2 * (omega_2 * omega_2) * sin(d_theta) +
               2 * g * m_2 * cos(theta_2) * sin(d_theta) +
               2 * g * m_1 * sin(theta_1) +
               gamma_1) /
              (-2 * L_1 * (m_1 + m_2 * sin(d_theta) * sin(d_theta))));

    f_2[2] = ((m_2 * L_2 * (omega_2 * omega_2) * sin(2 * d_theta) +
               2 * (m_1 + m_2) * L_1 * (omega_1 * omega_1) * sin(d_theta) +
               2 * g * (m_1 + m_2) * cos(theta_1) * sin(d_theta) +
               gamma_2) /
              (2 * L_2 * (m_1 + m_2 * sin(d_theta) * sin(d_theta))));

    std::vector<cpl::Vector> derivs = {f_1, f_2};
    
    return derivs;
}


int main(int argc, char* argv[]) {

    
    std::cout << " Nonlinear damped driven double pendulum\n"
              << " --------------------------------\n";

    std::string mode = argv[1];

    m_1 = atof(argv[2]);
    m_2 = atof(argv[3]);
    L_1 = atof(argv[4]);                // Length of pendulum L_2
    L_2 = atof(argv[5]);                // Length of pendulum L_1
    q_1 = atof(argv[6]);                // Damping coefficient q_1
    q_2 = atof(argv[7]);                // Damping coefficient q_2
    Omega_D_1 = atof(argv[8]);          // Driving frequency Omega_D_1
    Omega_D_2 = atof(argv[9]);          // Driving frequency Omega_D_2
    F_D_1 = atof(argv[10]);              // Driving amplitude F_D_1
    F_D_2 = atof(argv[11]);             // Driving amplitude F_D_2
    theta_1 = atof(argv[12]);           // Angle of deflection of pendulum (0_1)
    theta_2 = atof(argv[13]);           // Angle of deflection of pendulum (0_2)
    omega_1 = atof(argv[14]);           // Velocity of pendulum (ω_1)
    omega_2 = atof(argv[15]);           // Velocity of pendulum (ω_2)

    t_iteration = atof(argv[16]);             // Integration time t_iteration
    dt = atof(argv[17]);                // Stepsize
    accuracy = atof(argv[18]);          // Accuracy of simulation
    
    std::ofstream dataFile("..\\out\\pendulum_double.dat");
    std::ofstream coordinates("..\\out\\pendulum_double_coords.dat");

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
    dataFile << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2
                << '\t' << dt << '\t' << energy_1() << '\t' << energy_2() << '\t' << 0 << '\n';
    
    std::cout << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2
                << '\t' << dt << '\t' << energy_1() << '\t' << energy_2() << '\t' << 0 << '\n';

    coordinates << t << '\t' << L_1 * sin(theta_1) << '\t' << - L_1 * cos(theta_1) << '\t'
                             << L_1 * sin(theta_1) + L_2 * sin(theta_2) << '\t' << - (L_1 * cos(theta_1) + L_2 * cos(theta_2)) << '\n';

    while (t < t_iteration) {
        if(mode == "runge") {
            cpl::RK4Step_double(x, y, dt, derivates);
        }

        else if(mode == "rkck") {
            cpl::RKCKStep_double(x, y, dt, derivates);
        }

        else if(mode == "euler") {
            cpl::EulerStep_double(x, y, dt, derivates);
        }

        else if(mode == "eulercromer") {
            cpl::EulerCromerStep_double(x, y, dt, derivates);
        }

        t = x[0];
        theta_1 = x[1], omega_1 = x[2];
        theta_2 = y[1], omega_2 = y[2];
        error = dt;

        while (theta_1 >= pi) theta_1 -= 2 * pi;
        while (theta_1 < -pi) theta_1 += 2 * pi;
        while (theta_2 >= pi) theta_2 -= 2 * pi;
        while (theta_2 < -pi) theta_2 += 2 * pi;
        
        // Runtime test
        duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

        // Returns in radian (0) and radian\s (ω)
        dataFile << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2
                    << '\t' << error << '\t' << energy_1() << '\t' << energy_2() << '\t' << duration.count() << '\n';
        
        std::cout << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2
                    << '\t' << error << '\t' << energy_1() << '\t' << energy_2() << '\t' << duration.count() << '\n';

        coordinates << t << '\t' << L_1 * sin(theta_1) << '\t' << - L_1 * cos(theta_1) << '\t'
                                 << L_1 * sin(theta_1) + L_2 * sin(theta_2) << '\t' << - (L_1 * cos(theta_1) + L_2 * cos(theta_2)) << '\n';
    }

    std::cout << " Output data to file pendulum_double.dat" << std::endl;
    dataFile.close();
}

