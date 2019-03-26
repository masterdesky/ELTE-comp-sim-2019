#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

// Simulation parameters
int N;                    // Number of particles
double rho;               // Density (number per unit volume)
double T;                 // Temperature
double L;                 // Will be computed from N and rho

double boltzmann = 1.38e-23; // Boltzmann constant

double **r, **v, **a;     // Positions, velocities, accelerations
std::string boundary;

double Energy_current;    // Instantenous total energy
double Virial;            // Sum (r_ij * F_ij) in Virial
bool potential = false;   // Add potential energy to total energy

// Declare some functions
void initPositions();
void initVelocities();
void rescaleVelocities();
double instantaneousTemperature();

// Variables to implement Verlet's neighbor list
double rCutOff;           // Cut-off on Lennard-Jones potential and force
double rMax;              // Maximum separation to include in pair list
int nPairs;               // Number of pairs currently in pair list
int **pairList;           // The list of pair indices (i,j)
double **drPair;          // Vector separations of each pair (i,j)
double *rSqdPair;         // Squared separation of each pair (i,j)
int updateInterval     ;  // Number of time steps between updates of pair list

// Declare functions to implement neighbor list
void computeSeparation(int, int, double[], double&);
void updatePairList();
void updatePairSeparations();

void initialize() {
    r = new double* [N];
    v = new double* [N];
    a = new double* [N];
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        a[i] = new double [3];
    }
    initPositions();
    initVelocities();

    // Allocate memory for neighbor list variables
    nPairs = N * (N - 1) / 2;
    pairList = new int* [nPairs];
    drPair = new double* [nPairs];
    for (int p = 0; p < nPairs; p++) {
        pairList[p] = new int [2];      // To store indices i and j
        drPair[p] = new double [3];     // To store components x,y,z
    }
    rSqdPair = new double [nPairs];
}

void computeSeparation (int i, int j, double dr[], double& rSqd) {

    // Find separation using closest image convention
    rSqd = 0;
    for (int d = 0; d < 3; d++) {
        dr[d] = r[i][d] - r[j][d];
        if(boundary == "periodic") {
            // Closest image convention
            if (dr[d] >= 0.5*L) {
                dr[d] -= L;
            }
            if (dr[d] < -0.5*L) {
                dr[d] += L;
            }
        }
        rSqd += dr[d]*dr[d];
    }
}

void updatePairList() {
    nPairs = 0;
    double dr[3];
    for (int i = 0; i < N-1; i++)               // All distinct pairs
        for (int j = i+1; j < N; j++) {         // of particles i,j
            double rSqd;
            computeSeparation(i, j, dr, rSqd);
            if (rSqd < rMax*rMax) {
                pairList[nPairs][0] = i;
                pairList[nPairs][1] = j;
                ++nPairs;
            }
        }
}

void updatePairSeparations() {
    double dr[3];
    for (int p = 0; p < nPairs; p++) {
        int i = pairList[p][0];
        int j = pairList[p][1];
        double rSqd;
        computeSeparation(i, j, dr, rSqd);
        for (int d = 0; d < 3; d++)
            drPair[p][d] = dr[d];
        rSqdPair[p] = rSqd;
    }
}

void computeAccelerations() {

    for (int i = 0; i < N; i++)          // Set all accelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;

    for (int p = 0; p < nPairs; p++) {
        int i = pairList[p][0];
        int j = pairList[p][1];
        if (rSqdPair[p] < rCutOff*rCutOff) {
            double r2Inv = 1 / rSqdPair[p];
            double r6Inv = r2Inv*r2Inv*r2Inv;
            double f = 24*r2Inv*r6Inv*(2*r6Inv - 1);
            for (int d = 0; d < 3; d++) {
                a[i][d] += f * drPair[p][d];
                a[j][d] -= f * drPair[p][d];

                if(potential) {
                    if(i < j) {
                        Virial += abs(drPair[p][d]) * f;
                    }
                }
            }

            // Step with energy (+ potential energy)
            if(potential) {
                Energy_current += 4 * (pow(rSqdPair[p], -6) - pow(rSqdPair[p], -3));
            }
        }
    }
}

void velocityVerlet(double dt) {
    // Assume accelerations have been computed
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;

            if(boundary == "periodic") {
                // Use periodic boundary conditions
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
    updatePairSeparations();

    potential = true;
    computeAccelerations();
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {

            if(boundary == "bounded") {
                // use bounded boundary conditions
                if (r[i][k] < 0 || r[i][k] >= L) {
                    v[i][k] *= -1;
                    a[i][k] *= -1;
                }
            }

            v[i][k] += 0.5 * a[i][k] * dt;

            // Step with energy
            Energy_current += v[i][k]*v[i][k] * 0.5;
        }
    }
}

void initPositions() {

    // compute side of cube from number of particles and number density
    L = pow(N / rho, 1.0/3);

    // find M large enough to fit N atoms on an fcc lattice
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    double a = L / M;           // lattice constant of conventional cell

    // 4 atomic positions in fcc unit cell 
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};

    int n = 0;                  // atoms placed so far
    for (int x = 0; x < M; x++)
        for (int y = 0; y < M; y++)
            for (int z = 0; z < M; z++)
                for (int k = 0; k < 4; k++)
                    if (n < N) {
                        r[n][0] = (x + xCell[k]) * a;
                        r[n][1] = (y + yCell[k]) * a;
                        r[n][2] = (z + zCell[k]) * a;
                        ++n;
                    }
}

double gasdev () {
     static bool available = false;
     static double gset;
     double fac, rsq, v1, v2;
     if (!available) {
          do {
               v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
               rsq = v1 * v1 + v2 * v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = sqrt(-2.0 * log(rsq) / rsq);
          gset = v1 * fac;
          available = true;
          return v2*fac;
     } else {
          available = false;
          return gset;
     }
}

void initVelocities() {

    // Gaussian with unit variance
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev();
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];

    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities();
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

double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}

int main(int argc, char* argv[]) {

    boundary = argv[1];             // Mode for boundary conditions
    int n = atoi(argv[2]);          // Number of steps
    N = atoi(argv[3]);              // Number of particles
    rho = atof(argv[4]);            // Density (number per unit volume)
    T = atof(argv[5]);              // Temperature
    rCutOff = atof(argv[6]);        // Cut-off on Lennard-Jones potential and force
    rMax = atof(argv[7]);           // Maximum separation to include in pair list
    updateInterval = atoi(argv[8]); // Number of time steps between updates of pair list

    initialize();
    updatePairList();
    updatePairSeparations();
    computeAccelerations();
    double dt = 0.01;

    ofstream file("..\\out\\md3.dat");
    for (int i = 0; i < n; i++) {

        Energy_current = 0;
        Virial = 0;
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

        // Current energy
        file << Energy_current << '\t';

        // Virial
        file << Virial << '\t';

        // Temperature
        file << instantaneousTemperature() << '\n';

        if (i % 200 == 0)
            rescaleVelocities();
        if (i % updateInterval == 0) {
            updatePairList();
            updatePairSeparations();
        }
    }
    file.close();
}