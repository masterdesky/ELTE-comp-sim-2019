#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

// simulation parameters
int N;                    // number of particles
double rho;               // density (number per unit volume)
double T;                 // temperature
double dt;                // integration time step

void getInput() {
    cout << "Molecular Dynamics of 3D Lennard-Jones Gas" << endl
         << "------------------------------------------" << endl;
    cout << "Enter desired number of particles N = ";
    cin >> N;
    cout << "Enter density (particles/unit vol) rho = ";
    cin >> rho;
    cout << "Enter desired temperature T = ";
    cin >> T;
    cout << "Enter desired integration time step dt = ";
    cin >> dt;
}

double L;                 // will be computed from N and rho

double *r, *v, *a;        // positions, velocities, accelerations

// declare functions repeated from md-3D.cpp
void initPositions();
void initVelocities();
void rescaleVelocities();
double instantaneousTemperature();

// variables to implement Verlet's neighbor list
double rCutOff = 2.5;     // cut-off on Lennard-Jones potential and force
double rMax = 3.3;        // maximum separation to include in pair list
int nPairs;               // number of pairs currently in pair list
int **pairList;           // the list of pair indices (i,j)
double *drPair;           // vector separations of each pair (i,j)
double *rSqdPair;         // squared separation of each pair (i,j)
int updateInterval = 10;  // number of time steps between updates of pair list

// declare functions repeated from md3.cpp
void initialize();
void computeSeparation(int, int, double[], double&);
void updatePairList();
void updatePairSeparations();
void computeAccelerations();
void velocityVerlet(double dt);

int step;                 // keeps track of integration step number
int displayInterval = 5;  // display molecules every so many steps

void makeMolecules();     // this function re-draws the molecules

void takeStep() {
    velocityVerlet(dt);
    ++step;
    if (step % 200 == 0)
        rescaleVelocities();
    if (step % updateInterval == 0) {
        updatePairList();
        updatePairSeparations();
    }
    if (step % displayInterval == 0) {
        makeMolecules();
        glutPostRedisplay();
    }
}

const double pi = 4 * std::atan(1.0);
double radius = 0.5;                      // radius of molecule
double minExtent[3], maxExtent[3];        // extent of system volume
int xWindowSize = 640, yWindowSize = 640; // window size in screen pixels
GLdouble aspectRatio;                     // window aspect ratio
GLdouble fovy, nearClip, farClip;         // variables for 3D projection
GLdouble eye[3], center[3], up[3];        // more projection variables
GLuint sphereID, configID;                // display list ID numbers 
int phi, theta;                           // to rotate system using arrow keys
int angle = 5;                            // rotation angle in degrees

void makeSphere(GLuint listID, double radius) {
    int nTheta = 9;                       // number of polar angle slices
    int nPhi = 18;                        // number of azimuthal angle slices
    glNewList(listID, GL_COMPILE);
        glutSolidSphere(radius, nPhi, nTheta);
    glEndList();
}

void makeMolecules() {
    glNewList(configID, GL_COMPILE);
    glColor3f(1.0, 0.0, 0.0);             // color the molecules red
    glPushMatrix();
    glRotated(phi, 0, 1, 0);              // rotate about y axis
    glPushMatrix();
    glRotated(theta, 1, 0, 0);            // rotate about x axis
    for (int i = 0; i < N; i++) {
        glPushMatrix();
        glTranslated(r[3*i+0] - L/2, r[3*i+1] - L/2, r[3*i+2] - L/2);
        glCallList(sphereID);
        glPopMatrix();
    }
    glColor3ub(255, 255, 255);            // white
    glutWireCube(L);                      // cubical system volume
    glPopMatrix();
    glPopMatrix();
    glEndList();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2], 
              center[0], center[1], center[2],
              up[0], up[1], up[2]);
    glCallList(configID);                 // draw molecules
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    aspectRatio = w / double(h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovy, aspectRatio, nearClip, farClip);
    glMatrixMode(GL_MODELVIEW);
}

void initView(double *minExtent, double *maxExtent) {

    // use a single light source to illuminate the scene
    GLfloat lightDiffuse[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat lightPosition[] = {0.5, 0.5, 1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);    
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);              // use single light number 0
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);

    // compute the distance scale of the system
    double difExtent[3];
    for (int i = 0; i < 3; i++)
        difExtent[i] = maxExtent[i] - minExtent[i];
    double dist = 0;
    for (int i = 0; i < 3; i++)
        dist += difExtent[i] * difExtent[i];
    dist = std::sqrt(dist);

    // locate the center of the system, camera position, and orientation
    for (int i = 0; i < 3; i++)
        center[i] = minExtent[i] + difExtent[i] / 2;
    eye[0] = center[0];  
    eye[1] = center[1];  
    eye[2] = center[2] + dist;        // along z axis is towards viewer
    up[0] = 0;
    up[1] = 1;                        // y axis is up
    up[2] = 0;

    // set up clipping planes, field of view angle in degrees in y direction
    nearClip = (dist - difExtent[2] / 2) / 2;
    farClip = 2 * (dist + difExtent[2] / 2);
    fovy = difExtent[1] / (dist - difExtent[2] / 2) / 2;
    fovy = 2 * std::atan(fovy) / pi * 180;
    fovy *= 1.2;
}

void special(int key, int x, int y) {
    switch(key) {
      case GLUT_KEY_LEFT:   phi = (phi - angle) % 360;      break;
      case GLUT_KEY_RIGHT:  phi = (phi + angle) % 360;      break;
      case GLUT_KEY_UP:     theta = (theta - angle) % 360;  break;
      case GLUT_KEY_DOWN:   theta = (theta + angle) % 360;  break;
      default: break;
    }
}

bool running = false;

void mouse(int button, int state, int x, int y) {
    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running) {
                glutIdleFunc(NULL);
                running = false;
            } else {
                glutIdleFunc(takeStep);
                running = true;
            }
        }
    }
}

void menu(int heat) {
    double scale = heat ? 1.2 : 1/1.2;
    T *= scale;
    rescaleVelocities();
}

int main(int argc, char *argv[]) {

    glutInit(&argc, argv);

    getInput();
    initialize();
    updatePairList();
    updatePairSeparations();
    computeAccelerations();
    dt = 0.01;

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(xWindowSize, yWindowSize);
    glutCreateWindow("Molecular Dynamics of Lennard-Jones Gas");

    for (int i = 0; i < 3; i++) {
        minExtent[i] = -L/2;
        maxExtent[i] = L/2;
    } 
    initView(minExtent, maxExtent);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutSpecialFunc(special);
    glutMouseFunc(mouse);
    glutCreateMenu(menu);
    glutAddMenuEntry("Heat", 1);
    glutAddMenuEntry("Cool", 0);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    sphereID = glGenLists(1);
    makeSphere(sphereID, radius);
    configID = glGenLists(1);
    makeMolecules();

    glutMainLoop();
}

void initialize() {
    r = new double [3*N];
    v = new double [3*N];
    a = new double [3*N];
    initPositions();
    initVelocities();

    // allocate memory for neighbor list variables
    nPairs = N * (N - 1) / 2;
    pairList = new int* [nPairs];
    for (int p = 0; p < nPairs; p++)
        pairList[p] = new int [2];      // to store indices i and j
    drPair = new double [3*nPairs];
    rSqdPair = new double [nPairs];
}

void computeSeparation (int i, int j, double dr[], double& rSqd) {

    // find separation using closest image convention
    rSqd = 0;
    for (int d = 0; d < 3; d++) {
        dr[d] = r[3*i+d] - r[3*j+d];
        if (dr[d] >= 0.5*L)
            dr[d] -= L;
        if (dr[d] < -0.5*L)
            dr[d] += L;
        rSqd += dr[d]*dr[d];
    }
}

void updatePairList() {
    nPairs = 0;
    double *dr = new double [3];
    for (int i = 0; i < N-1; i++)               // all distinct pairs
        for (int j = i+1; j < N; j++) {         // of particles i,j
            double rSqd;
            computeSeparation(i, j, dr, rSqd);
            if (rSqd < rMax*rMax) {
                pairList[nPairs][0] = i;
                pairList[nPairs][1] = j;
                ++nPairs;
            }
        }
    delete [] dr;
}

void updatePairSeparations() {
    double *dr = new double [3];
    for (int p = 0; p < nPairs; p++) {
        int i = pairList[p][0];
        int j = pairList[p][1];
        double rSqd;
        computeSeparation(i, j, dr, rSqd);
        for (int d = 0; d < 3; d++)
            drPair[3*p+d] = dr[d];
        rSqdPair[p] = rSqd;
    }
    delete [] dr;
}

void computeAccelerations() {

    for (int i = 0; i < 3*N; i++)          // set all accelerations to zero
        a[i] = 0;

    for (int p = 0; p < nPairs; p++) {
        int i = pairList[p][0];
        int j = pairList[p][1];
        if (rSqdPair[p] < rCutOff*rCutOff) {
            double r2Inv = 1 / rSqdPair[p];
            double r6Inv = r2Inv*r2Inv*r2Inv;
            double f = 24*r2Inv*r6Inv*(2*r6Inv - 1);
            for (int d = 0; d < 3; d++) {
                a[3*i+d] += f * drPair[3*p+d];
                a[3*j+d] -= f * drPair[3*p+d];
            }
        }
    }
}

void velocityVerlet(double dt) {

    // assume accelerations have been computed and update positions
    for (int i = 0; i < 3*N; i++)
        r[i] += v[i] * dt + 0.5 * a[i] * dt * dt;

    // use periodic boundary conditions
    for (int i = 0; i < 3*N; i++) {
        if (r[i] < 0)
            r[i] += L;
        if (r[i] >= L)
            r[i] -= L;
    }

    for (int i = 0; i < 3*N; i++)
        v[i] += 0.5 * a[i] * dt;     // first half of velocity update

    updatePairSeparations();         // atoms have moved
    computeAccelerations();          // with updated separations

    for (int i = 0; i < 3*N; i++)
        v[i] += 0.5 * a[i] * dt;     // second half of velocity update
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
                        r[3*n+0] = (x + xCell[k]) * a;
                        r[3*n+1] = (y + yCell[k]) * a;
                        r[3*n+2] = (z + zCell[k]) * a;
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
    for (int i = 0; i < 3 * N; i++)
        v[i] = gasdev();
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[3*n+i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[3*n+i] -= vCM[i];

    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities();
}

void rescaleVelocities() {
    double vSqdSum = 0;
    for (int i = 0; i < 3 * N; i++)
        vSqdSum += v[i] * v[i];
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int i = 0; i < 3 * N; i++)
        v[i] *= lambda;
}

double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < 3 * N; i++)
        sum += v[i] * v[i];
    return sum / (3 * (N - 1));
}

