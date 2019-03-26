#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

const int N;                // Number of particles
//double r[N][3];           // Positions
//double v[N][3];           // Velocities
//double a[N][3];           // Accelerations
double T;                   // Temperature

double **r, **v, **a;     // positions, velocities, accelerations

double L = 10;            // Linear size of cubical volume
double vMax = 0.1;        // Maximum initial velocity component

double boltzmann = 1.38e-23; // Boltzmann constant

double Energy_current;    // Instantenous total energy
double Virial = 0;        // Sum (r_ij * F_ij) in Virial
bool potential;           // Add potential energy to total energy

std::string boundary;

void initialize() {

    r = new double* [N];
    v = new double* [N];
    a = new double* [N];
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        a[i] = new double [3];
    }

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

                if(potential) {
                    if(i < j) {
                        Virial += rij[k] * f;
                    }
                }
            }

            // Step with energy (+ potential energy)
            if(potential) {
                Energy_current += 4 * (pow(rSqd, -6) - pow(rSqd, -3));
            }
        }
}

void velocityVerlet(double dt) {
    potential = false;
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
    potential = true;
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

            // Step with energy (+ kinetic energy)
            Energy_current += v[i][k]*v[i][k] * 0.5;
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
    N = atoi(argv[3]);              // Number of particles
    T = atoi(argv[4]);              // Temperature

    double Energy = 0;              // Total energy
    double Energy2 = 0;             // Square of total energy

    initialize();
    double dt = 0.01;
    ofstream file("..\\out\\md1.dat");
    for (int i = 0; i < n; i++) {

        Energy_current = 0;
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

        double T_instant = instantaneousTemperature();

        // E; current energy
        file << Energy_current << '\t';

        Energy += Energy_current;
        Energy2 += Energy_current*Energy_current;

        // Average energy
        file << Energy/(i+1) << '\t';

        // Oscillation of energy
        double dE2 = ((Energy2)/(i+1) - (Energy/(i+1))*(Energy/(i+1)));
        file << dE2 << '\t';

        // C_v; molar heat capacity
        double C_v = 1/(boltzmann * T_instant * T_instant) * dE2;
        file << C_v << '\t';

        // PV
        double PV = N * boltzmann * T_instant + 1/3 * Virial/(i+1);
        file << PV / pow(L, 3) << '\t';

        // Z
        double Z = PV / (N * boltzmann * T_instant);
        file << Z << '\t';

        // Temperature
        file << T_instant << '\n';

        if (i % 200 == 0)
            rescaleVelocities();
    }
    file.close();
}

