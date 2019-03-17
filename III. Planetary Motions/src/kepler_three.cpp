#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>

#include "vector.hpp"
#include "odeint.hpp"

const double pi = 4 * atan(1.0);
const double GMPlusm = 4 * pi * pi;         // Kepler's Third Law: G(M + m)/(4*pi^2) = 1 [AU^3/year^2]
const double G = 1.9838 * pow(10, -29);     // Gravitational constant [AU^3 * kg^-1 * year^-2]
const double c = 63197.8;                   // Speed of light [AU/year]

static double m_1;                          // Mass of first body [kg]
static double m_2;                          // Mass of second body [kg]
static double m_3;                          // Mass of the third (perturbing) body [kg]
static double m_temp_1;                     // First of the two big masses for second derivate calculations [kg]
static double m_temp_2;                     // Second of the two big masses for second derivate calculations [kg]
static double r;                            // Distance of the two bodies' aphelions [AU]
static double r_ap_1;                       // Aphelion distance of first body [AU]
static double r_ap_2;                       // Aphelion distance of second body [AU]
static double eccentricity_1;               // Eccentricity of first body
static double eccentricity_2;               // Eccentricity of second body
static double a_1;                          // Length of semi-major axis of first body [AU]
static double a_2;                          // Length of semi-major axis of second body [AU]
static double v0_1;                         // Initial velocity of first body (tangential along y-axis) [AU/year]
static double v0_2;                         // Initial velocity of second body (tangential along y-axis) [AU/year]

std::vector<double> CurrentCoordinates(4);

static double plotting_years;               // Number of calculated years [year]
static double dt;                           // Step size [year]
static double accuracy;                     // Adaptive accuracy of simulation
static double number_of_bodies;             // Number of small bodies in the simulation

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
    double t = x[0], r_x_Curr = x[1], r_y_Curr = x[2], v_x_Curr = x[3], v_y_Curr = x[4];

    double rSquared_1 = (r_x_Curr - CurrentCoordinates[0])*(r_x_Curr - CurrentCoordinates[0]) + (r_x_Curr - CurrentCoordinates[1])*(r_x_Curr - CurrentCoordinates[3]);
    double rSquared_2 = (r_x_Curr - CurrentCoordinates[2])*(r_x_Curr - CurrentCoordinates[2]) + (r_x_Curr - CurrentCoordinates[1])*(r_x_Curr - CurrentCoordinates[3]);

    double rCubed_1 = rSquared_1 * sqrt(rSquared_1);
    double rCubed_2 = rSquared_2 * sqrt(rSquared_2);

    cpl::Vector f(5);
    f[0] = 1;
    f[1] = v_x_Curr;
    f[2] = v_y_Curr;
    f[3] = - G * (m_temp_1 * CurrentCoordinates[0]/rCubed_1 - m_temp_2 * CurrentCoordinates[2]/rCubed_2);
    f[4] = - G * (m_temp_1 * CurrentCoordinates[1]/rCubed_1 - m_temp_2 * CurrentCoordinates[3]/rCubed_2);

    // Relativistic effects for Keplerian orbit, due to special relativity
    if(relat) {
        double gamma_x = sqrt(1 - pow(v_x_Curr, 2)/pow(c, 2));
        double gamma_y = sqrt(1 - pow(v_y_Curr, 2)/pow(c, 2));
        f[3] /= gamma_x;
        f[4] /= gamma_y;
    }
    if(switch_t_with_y) {
        //  use y as independent variable
        for(int i = 0; i < 5; i++) {
            f[i] /= v_y_Curr;
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
    number_of_bodies = atof(argv[12]);              // Number of small bodies in the simulation

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

    std::vector<cpl::Vector> SmallBodies(number_of_bodies, cpl::Vector(5));
    std::vector<double> SmallBodyMasses(number_of_bodies);

    // Read in data of small bodies into a vector
    std::ifstream inputFile("small_objects.dat");

    int number_of_bodies = 0;                   // Number of small bodies in the simulation
    // Check if input file exists
    if(inputFile.good()) {
        
        double current_number = 0;              // Storage for current read number from data file
        int current_column = 0;
        int current_row = 0;
        // Push items into a vector
        while (inputFile >> current_number) {
            std::cout << current_number << '\t';

            if(current_column == 0) {
                SmallBodyMasses[current_row] = current_number;
                SmallBodies[current_row][0] = 0;
            }

            else{
                SmallBodies[current_row][current_column] = current_number;
            }

            current_column += 1;

            if(current_column == 5) {
                current_column = 0;
                current_row++;
                std::cout << '\n';
            }
        }

        // Close the file.
        inputFile.close();
    }

    else {
        exit(-1);
    }

    // Print read in datapoints
    std::cout << "Read in data:\n"
              << " -------------------------------\n";
    for(int i = 0; i < number_of_bodies; i++) {
        for(int j = 0; j < 5; j++) {
            std::cout << SmallBodies[i][j] << '\t';
        }
        std::cout << '\n';
    }

    //  Changing variables
    std::ofstream dataFile;         // Datafile for outputs
    std::ofstream dataFile_small;   // Datafile for small objects' outputs
    cpl::Vector x_1;                // Storing orbit parameters of first body for every step
    cpl::Vector x_2;                // Storing orbit parameters of second body for every step
    cpl::Vector x_3;                // Temporary storage for coordinates of small objects
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
        //  Fixed step size datafile and variable container for small bodies
        dataFile_small.open("..\\out\\fixed_smalls.dat");

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
                m_temp_1 = m_2;
                m_temp_2 = 0;
                cpl::RK4Step(x_1, dt, derivates);

                m_temp_1 = m_1;
                m_temp_2 = 0;
                cpl::RK4Step(x_2, dt, derivates);

                m_temp_1 = m_1;
                m_temp_2 = m_2;
                CurrentCoordinates[0] = x_1[1], CurrentCoordinates[1] = x_1[2], CurrentCoordinates[2] = x_2[1], CurrentCoordinates[3] = x_2[2];
                for(int i = 0; i < number_of_bodies; i++) {
                    x_3 = SmallBodies[i];
                    double y_3 = x_3[2];
                    cpl::RK4Step(x_3, dt, derivates);

                    if(y_3 * x_3[2] < 0) {
                        interpolate_crossing(x_3, crossing);
                    }

                    for(int i = 0; i < 5; i++) {
                        dataFile_small << x_3[i] << '\t';
                    }
                }
            }
            else if(odeint=="rkck") {
                m_temp_1 = m_2;
                m_temp_2 = 0;
                cpl::RKCKStep(x_1, dt, derivates);

                m_temp_1 = m_1;
                m_temp_2 = 0;
                cpl::RKCKStep(x_2, dt, derivates);

                m_temp_1 = m_1;
                m_temp_2 = m_2;
                CurrentCoordinates[0] = x_1[1], CurrentCoordinates[1] = x_1[2], CurrentCoordinates[2] = x_2[1], CurrentCoordinates[3] = x_2[2];
                for(int i = 0; i < number_of_bodies; i++) {
                    x_3 = SmallBodies[i];
                    double y_3 = x_3[2];
                    cpl::RKCKStep(x_3, dt, derivates);

                    double steps_temp = steps;
                    steps++;
                    if(y_3 * x_3[2] < 0) {
                        interpolate_crossing(x_3, crossing);
                    }
                    
                    for(int i = 0; i < 5; i++) {
                        dataFile_small << x_3[i] << '\t';
                    }
                    steps = steps_temp;
                }
            }
            steps++;
            if(y_1 * x_1[2] < 0) {
                m_temp_1 = m_2;
                m_temp_2 = 0;
                interpolate_crossing(x_1, crossing);
            }
            if(y_2 * x_2[2] < 0) {
                m_temp_1 = m_1;
                m_temp_2 = 0;
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
        //  Adaptive step size datafile and variable container for small bodies
        dataFile.open("..\\out\\adaptive_smalls.dat");

        x_1 = x0_1;
        x_2 = x0_2;
        steps = crossing = 0;
        double dt_max_1 = 0, dt_min_1 = 100;
        double dt_max_2 = 0, dt_min_2 = 100;
        
        std::vector<double> dt_max_3(number_of_bodies), dt_min_3(number_of_bodies);

        for(int i = 0; i < number_of_bodies; i++) {
            dt_max_3[i] = 0;
            dt_min_3[i] = 100;
        }

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
                m_temp_1 = m_2;
                m_temp_2 = 0;
                cpl::adaptiveRK4Step(x_1, dt, accuracy, derivates);

                m_temp_1 = m_1;
                m_temp_2 = 0;
                cpl::adaptiveRK4Step(x_2, dt, accuracy, derivates);

                m_temp_1 = m_1;
                m_temp_2 = m_2;
                CurrentCoordinates[0] = x_1[1], CurrentCoordinates[1] = x_1[2], CurrentCoordinates[2] = x_2[1], CurrentCoordinates[3] = x_2[2];
                for(int i = 0; i < number_of_bodies; i++) {
                    x_3 = SmallBodies[i];
                    double t_save_3 = x_3[0];
                    double y_3 = x_3[2];
                    cpl::adaptiveRK4Step(x_3, dt, accuracy, derivates);

                    double step_size_3 = x_3[0] - t_save_3;
                    double steps_temp = steps;
                    steps++;
                    if(step_size_3 < dt_min_3[i]) {
                        dt_min_3[i] = step_size_3;
                    }
                    if(step_size_3 < dt_max_3[i]) {
                        dt_max_3[i] = step_size_3;
                    }

                    if(y_3 * x_3[2] < 0) {
                        interpolate_crossing(x_3, crossing);
                    }
                    steps = steps_temp;
                }
            }
            else if(odeint=="rkck") {
                m_temp_1 = m_2;
                m_temp_2 = 0;
                cpl::adaptiveRKCKStep(x_1, dt, accuracy, derivates);

                m_temp_1 = m_1;
                m_temp_2 = 0;
                cpl::adaptiveRKCKStep(x_2, dt, accuracy, derivates);

                m_temp_1 = m_1;
                m_temp_2 = m_2;
                CurrentCoordinates[0] = x_1[1], CurrentCoordinates[1] = x_1[2], CurrentCoordinates[2] = x_2[1], CurrentCoordinates[3] = x_2[2];
                for(int i = 0; i < number_of_bodies; i++) {
                    x_3 = SmallBodies[i];
                    double t_save_3 = x_3[0];
                    double y_3 = x_3[2];
                    cpl::adaptiveRKCKStep(x_3, dt, accuracy, derivates);

                    double step_size_3 = x_3[0] - t_save_3;
                    double steps_temp = steps;
                    steps++;
                    if(step_size_3 < dt_min_3[i]) {
                        dt_min_3[i] = step_size_3;
                    }
                    if(step_size_3 < dt_max_3[i]) {
                        dt_max_3[i] = step_size_3;
                    }

                    if(y_3 * x_3[2] < 0) {
                        interpolate_crossing(x_3, crossing);
                    }
                    steps = steps_temp;
                }
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
                m_temp_1 = m_2;
                m_temp_2 = 0;
                interpolate_crossing(x_1, crossing);
            }
            if(y_2 * x_2[2] < 0) {
                m_temp_1 = m_1;
                m_temp_2 = 0;
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

