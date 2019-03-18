#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

const int N = 64;         // number of particles
double r[N][3];           // positions
double v[N][3];           // velocities
double a[N][3];           // accelerations

double L = 10;            // linear size of cubical volume
double vMax = 0.1;        // maximum initial velocity component

void initialize() {

    // initialize positions
    int n = int(ceil(pow(N, 1.0/3)));  // number of atoms in each direction
    double a = L / n;                  // lattice spacing
    int p = 0;                         // particles placed so far
    for (int x = 0; x < n; x++) 
        for (int y = 0; y < n; y++) 
            for (int z = 0; z < n; z++) {
                    if (p < N) {
                        r[p][0] = (x + 0.5) * a;
                        r[p][1] = (y + 0.5) * a;
                        r[p][2] = (z + 0.5) * a;
                    }
                    ++p;
            }

    // initialize velocities
    for (int p = 0; p < N; p++)
        for (int i = 0; i < 3; i++)
            v[p][i] = vMax * (2 * rand() / double(RAND_MAX) - 1);
}

void computeAccelerations() {

    for (int i = 0; i < N; i++)          // set all accelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int i = 0; i < N-1; i++)        // loop over all distinct pairs i,j
        for (int j = i+1; j < N; j++) { 
            double rij[3];               // position of i relative to j
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];    
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            for (int k = 0; k < 3; k++) {
                 a[i][k] += rij[k] * f;
                 a[j][k] -= rij[k] * f;
            }
        }
}

void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}

double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}

int main() {
    initialize();
    double dt = 0.01;
    ofstream file("T.data");
    for (int i = 0; i < 1000; i++) {
        velocityVerlet(dt);
        file << instantaneousTemperature() << '\n';
    }
    file.close();
}

