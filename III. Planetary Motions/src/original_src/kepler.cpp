#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

#include "vector.hpp"
#include "odeint.hpp"
using namespace cpl;

const double pi = 4 * atan(1.0);
const double GmPlusM = 4 * pi * pi;

bool switch_t_with_y = false;    //  to interpolate to y = 0

//  Derivative vector for Newton's law of gravitation
Vector f(const Vector& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = r_x*r_x + r_y*r_y;
    double rCubed = rSquared * sqrt(rSquared);
    Vector f(5);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[3] = - GmPlusM * r_x / rCubed;
    f[4] = - GmPlusM * r_y / rCubed;
    if (switch_t_with_y) {
        //  use y as independent variable
        for (int i = 0; i < 5; i++)
            f[i] /= v_y;
    }
    return f;
}

//  Change independent variable from t to y and step back to y = 0
void interpolate_crossing(Vector x, int& crossing) {
    ++crossing;
    switch_t_with_y = true;
    RK4Step(x, -x[2], f);
    cout << " crossing " << crossing << "\t t = " << x[0]
         << "\t x = " << x[1] << endl;
    switch_t_with_y = false;
}

int main() {
    cout << " Kepler orbit comparing fixed and adaptive Runge-Kutta\n"
         << " -----------------------------------------------------\n"
         << " Enter aphelion distance in AU, and eccentricity: ";
    double r_ap, eccentricity, a, T, v0;
    cin >> r_ap >> eccentricity;
    a = r_ap / (1 + eccentricity);
    T = pow(a, 1.5);
    v0 = sqrt(GmPlusM * (2 / r_ap - 1 / a));
    cout << " Enter number of periods, step size, and adaptive accuracy: ";
    double periods, dt, accuracy;
    cin >> periods >> dt >> accuracy;
    Vector x0(5);
    x0[0] = 0;  x0[1] = r_ap;  x0[2] = 0;  x0[3] = 0;  x0[4] = v0;

    ofstream dataFile("fixed.data");
    Vector x = x0;
    int steps = 0, crossing = 0;
    cout << "\n Integrating with fixed step size" << endl;
    do {
        for (int i = 0; i < 5; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        double y = x[2];
        RK4Step(x, dt, f);
        ++steps;
        if (y * x[2] < 0)
            interpolate_crossing(x, crossing);
    } while (x[0] < periods * T);
    cout << " number of fixed size steps = " << steps << endl;
    cout << " data in file fixed.data" << endl;
    dataFile.close();

    dataFile.open("adaptive.data");
    x = x0;
    steps = crossing = 0;
    double dt_max = 0, dt_min = 100;
    cout << "\n Integrating with adaptive step size" << endl;
    do {
        for (int i = 0; i < 5; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        double t_save = x[0];
        double y = x[2];
        adaptiveRK4Step(x, dt, accuracy, f);
        double step_size = x[0] - t_save;
        ++steps;
        if (step_size < dt_min) dt_min = step_size;
        if (step_size > dt_max) dt_max = step_size;
        if (y * x[2] < 0)
            interpolate_crossing(x, crossing);
    } while (x[0] < periods * T);
    cout << " number of adaptive steps = " << steps << endl;
    cout << " step size: min = " << dt_min << "  max = " << dt_max << endl;
    cout << " data in file adaptive.data" << endl;
    dataFile.close();
}

