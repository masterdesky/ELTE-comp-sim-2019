#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>

#include "vector.hpp"
#include "odeint.hpp"

const double pi = 4 * atan(1.0);
const double G = 1.9838 * pow(10, -29);     // Gravitational constant [AU^3 * kg^-1 * year^-2]
const double c = 63197.8;                   // Speed of light [AU/year]

static double m_1;                          // Mass of the first body [kg]
static double m_2;                          // Mass of the second body [kg]
static double m;                            // Current orbiting mass for second derivate calculations [kg]
static double r_ap;                         // Aphelion distance [AU]
static double eccentricity;                 // Eccentricity
static double a;                            // Length of semi-major axis [AU]
static double v0;                           // Initial velocity (tangential along y-axis) [AU/year]

static double plotting_years;               // Number of calculated years [year]
static double dt;                           // Step size [year]
static double accuracy;                     // Adaptive accuracy of simulation
static double step_size;                    // Current adaptive stepsize

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
    f[3] = - G * m_2 * (m_1 + m_2)/m * r_x / rCubed;
    f[4] = - G * m_2 * (m_1 + m_2)/m * r_y / rCubed;

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
    
    std::string fixed_or_not = argv[1];         // Fixed or adaptive
    odeint = argv[2];                           // ODE integration method
    std::string relativity = argv[3];           // Relativistic effects
    m_1 = atof(argv[4]);                        // Mass of the first body [kg]
    m_2 = atof(argv[5]);                        // Mass of the second body [kg]
    r_ap = atof(argv[6]);                       // Aphelion distance [AU]
    eccentricity = atof(argv[7]);               // Eccentricity
    plotting_years = atof(argv[8]);             // Number of calculated years [year]
    dt = atof(argv[9]);                         // Step size [year]
    accuracy = atof(argv[10]);                  // Adaptive accuracy of simulation

    a = r_ap / (1 + eccentricity);              // Length of semi-major axis [AU]
    v0 = sqrt(G * (m_1 + m_2) * (2 / r_ap - 1 / a));    // Initial velocity (tangential along y-axis) [AU/year]

    if(relativity[0] == 'r') {
        relat = true;
    }

    //  Initial parameters
    //  x0[0]: time; x0[1]: x coordinate; x0[2]: y coordinate; x0[3]: x velocity; x0[4]: y velocity
    cpl::Vector x0(5);
    x0[0] = 0;  x0[1] = r_ap;  x0[2] = 0;  x0[3] = 0;  x0[4] = v0;

    //  Changing variables
    std::ofstream dataFile;     // Datafile for outputs
    cpl::Vector x;              // Storing orbit parameters for every step
    int steps, crossing;        // Stepsize, Interpolate crossing

    //
    //  FIXED STEP SIZE
    //
    if(fixed_or_not == "fixed") {
        // Time measurement init and starts
        auto start = std::chrono::steady_clock::now();
        std::chrono::microseconds duration;
        
        //  Fixed step size datafile and variable container 'x'
        dataFile.open("..\\out\\fixed.dat");
        x = x0;
        steps = 0, crossing = 0;

        //  Fixed step size
        std::cout << "\n Integrating with fixed step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 5; i++) {
                dataFile << x[i] << '\t';
            }

            dataFile << gravitational_potential(x, m_1) << '\t' << duration.count() << '\n';
        
            double y = x[2];
            if(odeint=="runge") {
                m = m_2;
                cpl::RK4Step(x, dt, derivates);
            }
            else if(odeint=="rkck") {
                m = m_2;
                cpl::RKCKStep(x, dt, derivates);
            }
            
            steps++;
            if(y * x[2] < 0) {
                m = m_2;
                interpolate_crossing(x, crossing);
            }

        } while (x[0] < plotting_years);

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
        x = x0;
        steps = crossing = 0;
        double dt_max = 0, dt_min = 100;

        //  Adaptive step size
        std::cout << "\n Integrating with adaptive step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 5; i++) {
                dataFile << x[i] << '\t';
            }

            dataFile << gravitational_potential(x, m_1) << '\t' << step_size << '\t' << duration.count() << '\n';
            double t_save = x[0];
            double y = x[2];
            if(odeint=="runge") {
                m = m_2;
                cpl::adaptiveRK4Step(x, dt, accuracy, derivates);
            }
            else if(odeint=="rkck") {
                m = m_2;
                cpl::adaptiveRKCKStep(x, dt, accuracy, derivates);
            }
            
            step_size = x[0] - t_save;
            steps++;
            if(step_size < dt_min) {
                dt_min = step_size;
            }
            if(step_size > dt_max) {
                dt_max = step_size;
            }
            if(y * x[2] < 0) {
                m = m_2;
                interpolate_crossing(x, crossing);
            }
        } while (x[0] < plotting_years);

        std::cout << " number of adaptive steps = " << steps << std::endl;
        std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;
        std::cout << " data in file adaptive.dat" << std::endl;
        dataFile.close();
    }
}

