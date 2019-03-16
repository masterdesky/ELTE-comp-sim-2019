#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>

#include "vector.hpp"
#include "odeint.hpp"

const double pi = 4 * atan(1.0);
const double GMPlusm = 4 * pi * pi;         // Kepler's Third Law: G(M + m)/(4*pi^2) = 1 [AU^3/year^2]
const double G = 1.9838 * pow(10, -29);     // Gravitational constant [AU^3 * kg^-1 * year^-2]
const double c = 63197.8;                   // Speed of light [AU/year]

static double m_1;                          // Mass of first body [kg]
static double m_2;                          // Mass of second body [kg]
static double m;                            // Current orbiting mass for second derivate calculations [kg]
static double r;                            // Distance of the two bodies' aphelions [AU]
static double r_ap_1;                       // Aphelion distance of first body [AU]
static double r_ap_2;                       // Aphelion distance of second body [AU]
static double eccentricity_1;               // Eccentricity of first body
static double eccentricity_2;               // Eccentricity of second body
static double a_1;                          // Length of semi-major axis of first body [AU]
static double a_2;                          // Length of semi-major axis of second body [AU]
static double v0_1;                         // Initial velocity of first body (tangential along y-axis) [AU/year]
static double v0_2;                         // Initial velocity of second body (tangential along y-axis) [AU/year]

static double plotting_years;               // Number of calculated years [year]
static double dt;                           // Step size [year]
static double accuracy;                     // Adaptive accuracy of simulation

std::string odeint;                         // ODE integration method
bool switch_t_with_y = false;               // To interpolate to y = 0
bool relat = false;                         // Enable/disable relativistic corrections


//  Calculate kinetic energy
auto kinetic_energy(const cpl::Vector& x, double& m) {
    // Dimension: kg * AU^2 / year^2
    // Convert to J, by changing AU^2/year^2 to m^2/s^2 (pow at the end)
    return 1/2 * m * (pow(x[3], 2) + pow(x[4], 2)) * pow(4740.57172, 2);
}

//  Calculate potential energy
double gravitational_potential(const cpl::Vector& x, double& M) {
    // Dimension: kg * AU^2 / year^2
    // Convert to J, by changing AU^2/year^2 to m^2/s^2 (pow at the end)
    return M * G / (sqrt(pow(x[1], 2) + pow(x[2], 2))) * pow(4740.57172, 2);
}


//  Derivative vector for Newton's law of gravitation
cpl::Vector derivates(const cpl::Vector& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = r_x*r_x + r_y*r_y;
    double rCubed = rSquared * sqrt(rSquared);

    cpl::Vector f(5);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[3] = - G * pow(m, 3)/pow((m_1 + m_2), 2) * r_x / rCubed;
    f[4] = - G * pow(m, 3)/pow((m_1 + m_2), 2) * r_y / rCubed;

    // Relativistic effects for Keplerian orbit, due to special relativity
    if(relat) {
        double gamma_x = sqrt(1 - pow(v_x, 2)/pow(c, 2));
        double gamma_y = sqrt(1 - pow(v_y, 2)/pow(c, 2));
        f[3] /= gamma_x;
        f[4] /= gamma_y;
    }
    if(switch_t_with_y) {
        //  use y as independent variable
        for(int i = 0; i < 5; i++) {
            f[i] /= v_y;
        }
    }
    return f;
}


//  Change independent variable from t to y and step back to y = 0
void interpolate_crossing(cpl::Vector x, int& crossing) {
    crossing++;
    switch_t_with_y = true;
    if(odeint=="runge") {
        cpl::RK4Step(x, -x[2], derivates);
    }
    else if(odeint=="rkck") {
        cpl::RKCKStep(x, -x[2], derivates);
    }
    std::cout << " crossing " << crossing << '\t' << " t = " << x[0] << '\t' << " x = " << x[1] << std::endl;
    switch_t_with_y = false;
}


int main(int argc, char* argv[]) {
    std::cout << " Kepler orbit comparing fixed and adaptive Runge-Kutta\n"
              << " -----------------------------------------------------\n";
    
    std::string fixed_or_not = argv[1];             // Fixed or adaptive
    odeint = argv[2];                               // ODE integration method
    std::string relativity = argv[3];               // Relativistic effects
    m_1 = atof(argv[4]);                            // Mass of first body [kg]
    m_2 = atof(argv[5]);                            // Mass of second body [kg]
    r = atof(argv[6]);                              // Distance of the two bodies' aphelions [AU]
    eccentricity_1 = atof(argv[7]);                 // Eccentricity of first body
    eccentricity_2 = atof(argv[8]);                 // Eccentricity of second body
    plotting_years = atof(argv[9]);                 // Number of calculated years [year]
    dt = atof(argv[10]);                            // Step size [year]
    accuracy = atof(argv[11]);                      // Adaptive accuracy of simulation

    r_ap_1 = m_2/(m_1+m_2) * r;                     // Aphelium distance of first body [AU]
    r_ap_2 = m_1/(m_1+m_2) * r;                     // Aphelium distance of second body [AU]
    a_1 = r_ap_1 / (1 + eccentricity_1);            // Length of semi-major axis of first body [AU]
    a_2 = r_ap_2 / (1 + eccentricity_2);            // Length of semi-major axis of second body [AU]
    v0_1 = sqrt((G * pow(m_2, 3)/pow((m_1 + m_2), 2)) * (2 / r_ap_1 - 1 / a_1));     // Initial velocity of first body (tangential along y-axis) [AU/year]
    v0_2 = sqrt((G * pow(m_1, 3)/pow((m_1 + m_2), 2)) * (2 / r_ap_2 - 1 / a_2));     // Initial velocity of second body (tangential along y-axis) [AU/year]

    if(relativity[0] == 'r') {
        relat = true;
    }

    //  Initial parameters
    //  x0_x[0]: time; x0_x[1]: x coordinate; x0_x[2]: y coordinate; x0_x[3]: x velocity; x0_x[4]: y velocity
    cpl::Vector x0_1(5), x0_2(5);
    x0_1[0] = 0;  x0_1[1] = -r_ap_1;  x0_1[2] = 0;  x0_1[3] = 0;  x0_1[4] = v0_1;
    x0_2[0] = 0;  x0_2[1] = r_ap_2;  x0_2[2] = 0;  x0_2[3] = 0;  x0_2[4] = -v0_2;

    //  Changing variables
    std::ofstream dataFile;         // Datafile for outputs
    cpl::Vector x_1;                // Storing orbit parameters of first body for every step
    cpl::Vector x_2;                // Storing orbit parameters of second body for every step
    int steps, crossing;            // Stepsize, Interpolate crossing


    //
    //  FIXED STEP SIZE
    //
    if(fixed_or_not == "fixed") {
        // Time measurement init and starts
        auto start = std::chrono::steady_clock::now();
        std::chrono::microseconds duration;

        //  Fixed step size datafile and variable container 'x'
        dataFile.open("..\\out\\fixed.dat");
        x_1 = x0_1;
        x_2 = x0_2;
        steps = 0, crossing = 0;

        //  Fixed step size
        std::cout << "\n Integrating with fixed step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 5; i++) {
                dataFile << x_1[i] << '\t';
            }
            for(int i = 0; i < 5; i++) {
                dataFile << x_2[i] << '\t';
            }
            dataFile << kinetic_energy(x_1, m_1) << '\t' << kinetic_energy(x_2, m_2) << '\t' << duration.count() << '\n';
            double y_1 = x_1[2];
            double y_2 = x_2[2];
            if(odeint=="runge") {
                m = m_2;
                cpl::RK4Step(x_1, dt, derivates);

                m = m_1;
                cpl::RK4Step(x_2, dt, derivates);
            }
            else if(odeint=="rkck") {
                m = m_2;
                cpl::RKCKStep(x_1, dt, derivates);

                m = m_1;
                cpl::RKCKStep(x_2, dt, derivates);

            }
            steps++;
            if(y_1 * x_1[2] < 0) {
                m = m_2;
                interpolate_crossing(x_1, crossing);
            }
            if(y_2 * x_2[2] < 0) {
                m = m_1;
                interpolate_crossing(x_2, crossing);
            }

        } while (x_1[0] < plotting_years);

        std::cout << " number of fixed size steps = " << steps << std::endl;
        std::cout << " data in file fixed.dat" << std::endl;
        dataFile.close();
    }


    //
    //  ADAPTIVE STEP SIZE
    //
    else if(fixed_or_not == "adaptive") {
        // Time measurement init and starts
        auto start = std::chrono::steady_clock::now();
        std::chrono::microseconds duration;

        //  Adaptive step size datafile and variable container 'x'
        dataFile.open("..\\out\\adaptive.dat");
        x_1 = x0_1;
        x_2 = x0_2;
        steps = crossing = 0;
        double dt_max_1 = 0, dt_min_1 = 100;
        double dt_max_2 = 0, dt_min_2 = 100;

        //  Adaptive step size
        std::cout << "\n Integrating with adaptive step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 5; i++) {
                dataFile << x_1[i] << '\t';
            }
            for(int i = 0; i < 5; i++) {
                dataFile << x_2[i] << '\t';
            }

            dataFile << kinetic_energy(x_1, m_1) << '\t' << kinetic_energy(x_2, m_2) << '\t' << duration.count() << '\n';
            double t_save_1 = x_1[0];
            double t_save_2 = x_2[0];
            double y_1 = x_1[2];
            double y_2 = x_2[2];
            if(odeint=="runge") {
                m = m_2;
                cpl::adaptiveRK4Step(x_1, dt, accuracy, derivates);

                m = m_1;
                cpl::adaptiveRK4Step(x_2, dt, accuracy, derivates);
            }
            else if(odeint=="rkck") {
                m = m_2;
                cpl::adaptiveRKCKStep(x_1, dt, accuracy, derivates);

                m = m_1;
                cpl::adaptiveRKCKStep(x_2, dt, accuracy, derivates);

            }
            
            double step_size_1 = x_1[0] - t_save_1;
            double step_size_2 = x_2[0] - t_save_2;
            steps++;
            if(step_size_1 < dt_min_1) {
                dt_min_1 = step_size_1;
            }
            if(step_size_2 < dt_min_2) {
                dt_min_2 = step_size_2;
            }

            if(step_size_1 > dt_max_1) {
                dt_max_1 = step_size_1;
            }
            if(step_size_2 > dt_max_2) {
                dt_max_2 = step_size_2;
            }

            if(y_1 * x_1[2] < 0) {
                m = m_2;
                interpolate_crossing(x_1, crossing);
            }
            if(y_2 * x_2[2] < 0) {
                m = m_1;
                interpolate_crossing(x_2, crossing);
            }
        } while (x_1[0] < plotting_years);

        std::cout << " number of adaptive steps = " << steps << std::endl;
        std::cout << " step size: min_1 = " << dt_min_1 << "  max_1 = " << dt_max_1 << std::endl;
        std::cout << " step size: min_2 = " << dt_min_2 << "  max_2 = " << dt_max_2 << std::endl;
        std::cout << " data in file adaptive.dat" << std::endl;
        dataFile.close();
    }
}

