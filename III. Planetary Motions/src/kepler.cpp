#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "vector.hpp"
#include "odeint.hpp"

const double pi = 4 * atan(1.0);
const double GMPlusm = 4 * pi * pi;         // Kepler's Third Law: G(M + m)/(4*pi^2) = 1 [AU^3/day^2]
const double G = 1.488 * pow(10, -34);      // Gravitational constant [AU^3 * kg^-1 * day^-2]
const double M_Sun = 1.989 * pow(10, 30);   // Mass of Sun [kg]
const double c = 173.145;                   // Speed of light [AU/day]

static double r_ap;                         // Aphelion distance [AU]
static double eccentricity;                 // Eccentricity
static double a;                            // Length of semi-major axis [AU]
static double v0;                           // Initial velocity (tangential along y-axis) [AU/day]

static double plotting_years;               // Number of calculated years [year]
static double dt;                           // Step size [year]
static double accuracy;                     // Adaptive accuracy of simulation

bool switch_t_with_y = false;               // To interpolate to y = 0
bool relat = false;                         // Enable/disable relativistic corrections

//  Calculate potential energy
double gravitational_potential(const cpl::Vector& x) {
    return M_Sun * G / (sqrt(pow(x[1], 2) + pow(x[2], 2)));
}


//  Derivative vector for Newton's law of gravitation
cpl::Vector derivates(const cpl::Vector& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = pow(r_x, 2) + pow(r_y, 2);
    double rCubed = pow(rSquared, 3/2);
    double gamma_x = 1 - pow(v_x, 2)/pow(c, 2);
    double gamma_y = 1 - pow(v_y, 2)/pow(c, 2);
    cpl::Vector f(5);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[3] = - GMPlusm * r_x / rCubed;
    f[4] = - GMPlusm * r_y / rCubed;

    if(relat) {
        f[3] /= gamma_x;               // Relativistic effects for Keplerian
        f[4] /= gamma_y;               // orbit, due to special relativity
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
    cpl::RK4Step(x, -x[2], derivates);
    std::cout << " crossing " << crossing << '\t' << " t = " << x[0] << '\t' << " x = " << x[1] << std::endl;
    switch_t_with_y = false;
}


int main(int argc, char* argv[]) {
    std::cout << " Kepler orbit comparing fixed and adaptive Runge-Kutta\n"
              << " -----------------------------------------------------\n";
    
    std::string fixed_or_not = argv[1];         // Fixed or adaptive
    std::string relativity = argv[2];           // Relativistic effects
    r_ap = atof(argv[3]);                       // Aphelion distance [AU]
    eccentricity = atof(argv[4]);               // Eccentricity
    plotting_years = atof(argv[5]);             // Number of calculated years [year]
    dt = atof(argv[6]);                         // Step size [year]
    accuracy = atof(argv[7]);                   // Adaptive accuracy of simulation

    a = r_ap / (1 + eccentricity);              // Length of semi-major axis [AU]
    v0 = sqrt(GMPlusm * (2 / r_ap - 1 / a));    // Initial velocity (tangential along y-axis) [AU/day]

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
        //  Fixed step size datafile and variable container 'x'
        dataFile.open("..\\out\\fixed.dat");
        x = x0;
        steps = 0, crossing = 0;

        //  Fixed step size
        std::cout << "\n Integrating with fixed step size" << std::endl;
        do {
            for(int i = 0; i < 5; i++) {
                dataFile << x[i] << '\t';
            }
            dataFile << gravitational_potential(x) << '\n';
            double y = x[2];
            cpl::RK4Step(x, dt, derivates);
            steps++;
            if(y * x[2] < 0) {
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
        //  Adaptive step size datafile and variable container 'x'
        dataFile.open("..\\out\\adaptive.dat");
        x = x0;
        steps = crossing = 0;
        double dt_max = 0, dt_min = 100;

        //  Adaptive step size
        std::cout << "\n Integrating with adaptive step size" << std::endl;
        do {
            for(int i = 0; i < 5; i++) {
                dataFile << x[i] << '\t';
            }

            dataFile << gravitational_potential(x) << '\n';
            double t_save = x[0];
            double y = x[2];
            cpl::adaptiveRK4Step(x, dt, accuracy, derivates);
            double step_size = x[0] - t_save;
            steps++;
            if(step_size < dt_min) {
                dt_min = step_size;
            }
            if(step_size > dt_max) {
                dt_max = step_size;
            }
            if(y * x[2] < 0) {
                interpolate_crossing(x, crossing);
            }
        } while (x[0] < plotting_years);

        std::cout << " number of adaptive steps = " << steps << std::endl;
        std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;
        std::cout << " data in file adaptive.dat" << std::endl;
        dataFile.close();
    }
}

