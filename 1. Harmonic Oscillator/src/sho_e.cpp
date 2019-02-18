#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

double omega;          // the natural frequency
double x, v;           // position and velocity at time t
int periods;           // number of periods to integrate
int stepsPerPeriod;    // number of time steps dt per period
std::string fileName;       // name of output file

void Euler(double dt);     // takes an Euler step
double energy();           // computes the energy

void Euler (double dt) {
    double a = - omega * omega * x;
    double v_temp = v;
    v += a * dt;
    x += v_temp * dt;
}

double energy ( ) {
    return 0.5 * (v * v + omega * omega * x * x);
}

int main(int argc, char* argv[]) {
    omega = atof(argv[1]);
    x = atof(argv[2]);
    v = atof(argv[3]);
    periods = atoi(argv[4]);
    stepsPerPeriod = atoi(argv[5]);
    fileName = "../out/sho_e.dat";

    std::ofstream file(fileName.c_str());
    if (!file) {
        std::cerr << "Cannot open " << fileName << "\nExiting ...\n";
        return 1;
    }
    const double pi = 4 * atan(1.0);
    double T = 2 * pi / omega;
    double dt = T / stepsPerPeriod;
    double t = 0;

    // Time measurement init and starts
    auto start = std::chrono::steady_clock::now();
    std::chrono::microseconds duration;

    file << t << '\t' << x << '\t' << v << '\t' << energy() << '\t' << duration.count() << '\n';
    for (int p = 1; p <= periods; p++) {
        for (int s = 0; s < stepsPerPeriod; s++) {
            Euler(dt);
            t += dt;
            duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
            file << t << '\t' << x << '\t' << v << '\t' << energy() << '\t' << duration.count() << '\n';
        }
        std::cout << "Period = " << p << "\tt = " << t
             << "\tx = " << x << "\tv = " << v
             << "\tenergy = " << energy()
             << "\truntieme = " << duration.count() << std::endl;
    }
    file.close();
}