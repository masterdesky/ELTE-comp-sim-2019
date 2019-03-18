#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include "vector.hpp"
using namespace cpl;

int N = 16;                      // number of molecules
double L = 6;                    // linear size of square region
double kT = 1;                   // initial kinetic energy/molecule
int skip = 0;                    // steps to skip to speed graphics

void getInput() {
    cout << "Molecular Dynamics simulation of 2D Lennard-Jones system\n";
    cout << "Enter number of molecules N = ";
    cin >> N;
    cout << "Enter size of square region L = ";
    cin >> L;
    cout << "Enter kinetic energy/molecule in units of kT: ";
    cin >> kT;
    cout << "Enter steps to skip to speed graphics: ";
    cin >> skip;
}

Vector x, y, vx, vy;        // position and velocity components
Vector ax, ay;              // acceleration components

void computeAccelerations() { 

    for (int i = 1; i <= N; i++)
        ax[i] = ay[i] = 0;

    for (int i = 1; i < N; i++)
    for (int j = i + 1; j <= N; j++) {
        double dx = x[i] - x[j];
        double dy = y[i] - y[j];
        // use closest periodic image
        if (abs(dx) > 0.5 * L)
            dx *= 1 - L / abs(dx);
        if (abs(dy) > 0.5 * L)
            dy *= 1 - L / abs(dy);
        double dr = sqrt(dx * dx + dy * dy);
        double f = 48 * pow(dr, -13.0) - 24 * pow(dr, -7.0);
        ax[i] += f * dx / dr;
        ay[i] += f * dy / dr;
        ax[j] -= f * dx / dr;
        ay[j] -= f * dy / dr;
    }
}

double t = 0;                   // time
double dt = 0.01;               // integration time step
int step = 0;                   // step number
double T;                       // temperature
double Tsum;                    // to compute average T
int step0;                      // starting step for computing average

int nBins = 50;                 // number of velocity bins
Vector vBins;                   // for Maxwell-Boltzmann distribution
double vMax = 4;                // maximum velocity to bin
double dv = vMax / nBins;       // bin size

void resetHistogram() {
    for (int bin = 0; bin <= nBins; bin++)
        vBins[bin] = 0;
    Tsum = 0;
    step0 = step;
}

void initialize() {

    // allocate storage
    x = Vector(N+1);
    y = Vector(N+1);
    vx = Vector(N+1);
    vy = Vector(N+1);
    ax = Vector(N+1);
    ay = Vector(N+1);
    vBins = Vector(nBins+1);

    // initialize positions on a square lattice
    int m = int(ceil(sqrt(double(N))));  // molecules in each direction
    double d = L / m;                    // lattice spacing
    int n = 0;                           // molecules placed so far
    for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
        if (n < N) {
            ++n;
            x[n] = (i + 0.5) * d;
            y[n] = (j + 0.5) * d;
        }
    }

    // initialize velocities with random directions
    double pi = 4 * atan(1.0);
    double v = sqrt(2 * kT);
    for (int n = 1; n <= N; n++) {
        double theta = 2 * pi * rand() / double(RAND_MAX);
        vx[n] = v * cos(theta);
        vy[n] = v * sin(theta);
    }
    
    computeAccelerations();
    resetHistogram();
    T = kT;
}

void timeStep () {

    t += dt;
    ++step;
    double K = 0;
    for (int i = 1; i <= N; i++) {

        // integrate using velocity Verlet algorithm
        x[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;
        y[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;

        // periodic boundary conditions
        if (x[i] < 0) x[i] += L;
        if (x[i] > L) x[i] -= L;
        if (y[i] < 0) y[i] += L;
        if (y[i] > L) y[i] -= L;

        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;

        computeAccelerations();

        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;

        double v = sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
        K += 0.5 * v * v;
        int bin = (int) (nBins * v / vMax);
        if (bin > 0 && bin <= nBins)
            ++vBins[bin];
    }
    Tsum += K / N;
    T = Tsum / (step - step0);
}

void scaleVelocities(int heat) {
    double scale = heat ? 1.2 : 1/1.2;
    for (int n = 1; n <= N; n++) {
        vx[n] *= scale;
        vy[n] *= scale;        
    }
    resetHistogram();
}

int mainWindow, molsWindow, histWindow;
int margin = 10;

void takeStep() {
    for (int i = 0; i < skip; i++)
        timeStep();
    timeStep();
    glutSetWindow(molsWindow);
    glutPostRedisplay();
    if (step % 10 == 0) {
        glutSetWindow(histWindow);
        glutPostRedisplay();
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
}

void drawText(const string& str, double x, double y) {
    glRasterPos2d(x, y);
    int len = str.find('\0');
    for (int i = 0; i < len; i++)
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]); 
}

void displayMols() {
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3ub(0, 0, 255);
    double pi = 4 * atan(1.0);    
    for (int n = 1; n <= N; n++) {
        glBegin(GL_TRIANGLE_FAN);
            glVertex2d(x[n], y[n]);
            double phi = 2 * pi / 24;
            for (int j = 0; j < 25; j++)
                glVertex2d(x[n] + 0.5 * cos(phi*j), 
                           y[n] + 0.5 * sin(phi*j));
        glEnd();    
    }
    glColor3ub(255, 255, 255);
    ostringstream os;
    os << "Step No: " << step << "   Time t = " << t << ends;
    drawText(os.str(), L / 50, L / 50);
    os.seekp(0);
    glutSwapBuffers();
}

double MaxwellBoltzmann(double v, double T) {
    double Pmax = sqrt(T) * exp(-0.5);
    return v * exp( -0.5 * v * v / T) / Pmax;
}

void displayHist() {
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3ub(0, 0, 0);
    ostringstream os;
    os << "Temperature T = " << T << ends;
    drawText(os.str(), vMax / 2, 0.95);
    os.seekp(0);
    os << "Press left button to heat/cool" << ends;
    drawText(os.str(), vMax / 2, 0.9);
    os.seekp(0);

    double Pmax = 0;
    for (int b = 1; b <= nBins; b++)
        if (Pmax < vBins[b])
           Pmax = vBins[b];
    if (Pmax == 0) {
        glutSwapBuffers();
        return;
    }
 
    // draw velocity histogram
    glColor3ub(0, 255, 255);
    double dv = vMax / nBins; 
    for (int b = 1; b <= nBins; b++) {
        double v = b * dv;
        double P = vBins[b] / Pmax;
        glBegin(GL_POLYGON);
            glVertex2d(v - dv/3, 0);
            glVertex2d(v - dv/3, P);
            glVertex2d(v + dv/3, P);
            glVertex2d(v + dv/3, 0);
        glEnd();
    }
 
    // compare with Maxwell Boltzmann
    glColor3ub(255, 0, 255);
    glBegin(GL_LINE_STRIP);
        glVertex2d(0, 0);
        for (int b = 1; b < nBins; b++) {
            double v = b * dv;
            double P = MaxwellBoltzmann(v, T);
            glVertex2d(v, P);
        }
    glEnd();

    glColor3ub(0, 0, 0);
    os << "v_max = " << vMax << ends;
    drawText(os.str(), 0.8 * vMax, 0.02);
    os.seekp(0);

    glutSwapBuffers();
}

void reshape(int w, int h) {

    glutSetWindow(molsWindow);
    int width = h - 2*margin;
    glutReshapeWindow(width, width);
    glViewport(0, 0, width, width);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, L, 0, L);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glutSetWindow(histWindow);
    glutPositionWindow(h + margin, margin);
    width = w - h - 2*margin;
    if (width < h - 2*margin)
        width = h - 2*margin;
    glutReshapeWindow(width, h - 2*margin);
    glViewport(0, 0, width, width);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, vMax, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

bool running = false;

void molsMouse(int button, int state, int x, int y) {
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

void makeWindows() {

    // main window
    glutInitWindowSize(800, 400);
    glutInitWindowPosition(100, 100);
    mainWindow = glutCreateWindow("Molecular Dynamics of 2D "
                                  " Lennard-Jones System");
    glShadeModel(GL_FLAT);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    // molecules window
    molsWindow = glutCreateSubWindow(mainWindow, margin, margin, 
                                     400 - 2*margin, 400 - 2*margin);
    glClearColor(1.0, 0.0, 0.0, 0.0);
    glutDisplayFunc(displayMols);
    glutMouseFunc(molsMouse);

    // histogram window
    histWindow = glutCreateSubWindow(mainWindow, 400 + margin, margin,
                                     400 - 2*margin, 400 - 2*margin);
    glClearColor(0.0, 1.0, 0.0, 0.0);
    glutDisplayFunc(displayHist);
    glutCreateMenu(scaleVelocities);
    glutAddMenuEntry("Heat", 1);
    glutAddMenuEntry("Cool", 0);
    glutAttachMenu(GLUT_LEFT_BUTTON);
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    if (argc == 1)
        getInput();
    initialize();   
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    makeWindows();
    glutMainLoop();
}

