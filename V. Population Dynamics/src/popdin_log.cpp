#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>

#include "vector.hpp"
#include "odeint.hpp"

static double n_0;                          // Starting number of population
static double k;                            // Maximal number of population
static double w_out;                        // Dying rate of logistic simulation
static double w_in;                         // Birth rate of logistic simulation

static double sim_time;                     // Simulated time [units]
static double dt;                           // Step size [units]
static double accuracy;                     // Adaptive accuracy of simulation
static double step_size;                    // Current adaptive stepsize

std::string odeint;                         // ODE integration method

//  Derivative vector for Newton's law of gravitation
cpl::Vector derivates(const cpl::Vector& x) {

    double t = x[0], n = x[1];

    cpl::Vector f(2);
    f[0] = 1;
    f[1] = (-w_out + w_in) * n * (1 - n/k);

    return f;
}

int main(int argc, char* argv[]) {
    std::cout << " Population dynamics simulations with logistic model\n"
              << " ---------------------------------------------------\n";
    
    std::string fixed_or_not = argv[1];         // Fixed or adaptive
    odeint = argv[2];                           // ODE integration method

    n_0 = atof(argv[3]);                        // Starting number of population
    k = atof(argv[4]);                          // Maximal number of population
    w_out = atof(argv[5]);                      // Dying rate of logistic simulation
    w_in = atof(argv[6]);                       // Birth rate of logistic simulation

    sim_time = atof(argv[7]);                   // Simulated time [units]
    dt = atof(argv[8]);                         // Step size [units]
    accuracy = atof(argv[9]);                   // Adaptive accuracy of simulation

    //  Initial parameters
    //  x0[0]: time; x0[1]: starting population;
    cpl::Vector x0(2);
    x0[0] = 0;  x0[1] = n_0;

    //  Changing variables
    std::ofstream dataFile;     // Datafile for outputs
    cpl::Vector x;              // Storing orbit parameters for every step
    int steps;                  // Stepsize

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
        steps = 0;

        //  Fixed step size
        std::cout << "\n Integrating with fixed step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 2; i++) {
                dataFile << x[i] << '\t';
            }

            dataFile << duration.count() << '\n';
        
            if(odeint=="runge") {
                cpl::RK4Step(x, dt, derivates);
            }
            else if(odeint=="rkck") {
                cpl::RKCKStep(x, dt, derivates);
            }
            else if(odeint=="euler") {
                cpl::EulerStep(x, dt, derivates);
            }
            else if(odeint=="eulercromer") {
                cpl::EulerCromerStep(x, dt, derivates);
            }
            
            steps++;

        } while (x[0] < sim_time);

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
        steps = 0;
        double dt_max = 0, dt_min = 100;

        //  Adaptive step size
        std::cout << "\n Integrating with adaptive step size" << std::endl;
        do {
            // Runtime test
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);

            for(int i = 0; i < 2; i++) {
                dataFile << x[i] << '\t';
            }

            dataFile << step_size << '\t' << duration.count() << '\n';
            double t_save = x[0];

            if(odeint=="runge") {
                cpl::adaptiveRK4Step(x, dt, accuracy, derivates);
            }
            else if(odeint=="rkck") {
                cpl::adaptiveRKCKStep(x, dt, accuracy, derivates);
            }
            // No adaptove methods are available
            // Falling back to fixed stepsize!
            else if(odeint=="euler") {
                cpl::EulerStep(x, dt, derivates);
            }
            else if(odeint=="eulercromer") {
                cpl::EulerCromerStep(x, dt, derivates);
            }
            
            step_size = x[0] - t_save;
            steps++;
            if(step_size < dt_min) {
                dt_min = step_size;
            }
            if(step_size > dt_max) {
                dt_max = step_size;
            }

        } while (x[0] < sim_time);

        std::cout << " number of adaptive steps = " << steps << std::endl;
        std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;
        std::cout << " data in file adaptive.dat" << std::endl;
        dataFile.close();
    }
}