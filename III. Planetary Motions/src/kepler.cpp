#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "vector.hpp"
#include "odeint.hpp"

const double pi = 4 * atan(1.0);
const double GmPlusM = 4 * pi * pi;

static double r_ap;
static double eccentricity;
static double a;
static double T;
static double v0;
static double periods;
static double dt;
static double accuracy;

bool switch_t_with_y = false;    //  to interpolate to y = 0

//  Derivative vector for Newton's law of gravitation
cpl::Vector derivates(const cpl::Vector& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = r_x*r_x + r_y*r_y;
    double rCubed = rSquared * sqrt(rSquared);
    cpl::Vector f(5);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[3] = - GmPlusM * r_x / rCubed;
    f[4] = - GmPlusM * r_y / rCubed;
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
    
    r_ap = atof(argv[1]);                       // Aphelion distance in AU
    eccentricity = atof(argv[2]);               // Eccentricity
    periods = atof(argv[3]);                    // Number of periods
    dt = atof(argv[4]);                         // Step size
    accuracy = atof(argv[5]);                   // Adaptive accuracy of simulation

    a = r_ap / (1 + eccentricity);              // Length of semi-major axis
    T = pow(a, 1.5);                            // Period T
    v0 = sqrt(GmPlusM * (2 / r_ap - 1 / a));    // Initial velocity (tangential along y-axis)

    // Initial parameters
    // x0[0]: time; x0[1]: x coordinate; x0[2]: y coordinate; x0[3]: x velocity; x0[4]: y velocity
    cpl::Vector x0(5);
    x0[0] = 0;  x0[1] = r_ap;  x0[2] = 0;  x0[3] = 0;  x0[4] = v0;

    // Changing variables
    std::ofstream dataFile;     // Datafile for outputs
    cpl::Vector x;              // Storing orbit parameters for every step
    int steps, crossing;        // Stepsize, Interpolate crossing

    //
    // FIXED STEP SIZE
    //
    // Fixed step size datafile and variable container 'x'
    dataFile.open("..\\out\\fixed.dat");
    x = x0;
    steps = 0, crossing = 0;

    // Fixed step size
    std::cout << "\n Integrating with fixed step size" << std::endl;
    do {
        for(int i = 0; i < 5; i++) {
            dataFile << x[i] << '\t';
        }
        dataFile << '\n';
        double y = x[2];
        cpl::RK4Step(x, dt, derivates);
        steps++;
        if(y * x[2] < 0) {
            interpolate_crossing(x, crossing);
        }

    } while (x[0] < periods * T);
    std::cout << " number of fixed size steps = " << steps << std::endl;
    std::cout << " data in file fixed.dat" << std::endl;
    dataFile.close();


    //
    // ADAPTIVE STEP SIZE
    //
    // Adaptive step size datafile and variable container 'x'
    dataFile.open("..\\out\\adaptive.dat");
    x = x0;
    steps = crossing = 0;
    double dt_max = 0, dt_min = 100;

    // Adaptive step size
    std::cout << "\n Integrating with adaptive step size" << std::endl;
    do {
        for(int i = 0; i < 5; i++) {
            dataFile << x[i] << '\t';
        }

        dataFile << '\n';
        double t_save = x[0];
        double y = x[2];
        cpl::adaptiveRK4Step(x, dt, accuracy, derivates);
        double step_size = x[0] - t_save;
        ++steps;
        if(step_size < dt_min) {
            dt_min = step_size;
        }
        if(step_size > dt_max) {
            dt_max = step_size;
        }
        if(y * x[2] < 0) {
            interpolate_crossing(x, crossing);
        }
    } while (x[0] < periods * T);
    std::cout << " number of adaptive steps = " << steps << std::endl;
    std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;
    std::cout << " data in file adaptive.dat" << std::endl;
    dataFile.close();
}

