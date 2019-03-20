#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

const int N = 64;         // Number of particles
double r[N][3];           // Positions
double v[N][3];           // Velocities
double a[N][3];           // Accelerations
double T;                 // Temperature

double L = 10;            // Linear size of cubical volume
double vMax = 0.1;        // Maximum initial velocity component

std::string boundary;

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

void rescaleVelocities() {
    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
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
                if(boundary == "periodic") {
                    // closest image convention   
                    if (abs(rij[k]) > 0.5 * L) {
                        if (rij[k] > 0)
                            rij[k] -= L;
                        else
                            rij[k] += L;
                    }
                }
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
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;

            if(boundary == "periodic") {
                // use periodic boundary conditions
                if (r[i][k] < 0) {
                    r[i][k] += L;
                }
                if (r[i][k] >= L) {
                    r[i][k] -= L;
                }
            }
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    }
    computeAccelerations();
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            if(boundary == "bounded") {
                // use periodic boundary conditions
                if (r[i][k] < 0 || r[i][k] >= L) {
                    v[i][k] *= -1;
                    a[i][k] *= -1;
                }
            }
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    }
}


double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            sum += v[i][k] * v[i][k];
        }
    }
    return sum / (3 * (N - 1));
}

int main(int argc, char* argv[]) {
    
    boundary = argv[1];             // Mode for boundary conditions
    int n = atoi(argv[2]);          // Number of steps
    T = atoi(argv[3]);              // Temperature

    initialize();
    double dt = 0.01;
    ofstream file("..\\out\\md1.dat");
    for (int i = 0; i < n; i++) {
        velocityVerlet(dt);
        
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < 3; k++) {
                file << r[j][k] << '\t';
            }
            for(int k = 0; k < 3; k++) {
                file << v[j][k] << '\t';
            }
            for(int k = 0; k < 3; k++) {
                file << a[j][k] << '\t';
            }
        }
        file << instantaneousTemperature() << '\n';
        if (i % 200 == 0)
            rescaleVelocities();
    }
    file.close();
}

