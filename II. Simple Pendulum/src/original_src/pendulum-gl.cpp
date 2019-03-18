#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#ifdef __APPLE__             // Mac OS X uses different header
#  include <GLUT/glut.h>
#else                        // Unix and Windows
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#include "vector.hpp"        // vectors with components of type double
#include "odeint.hpp"        // ODE integration routines, Runge-Kutta ...
using namespace cpl;

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

double L = 1.0;              // length of pendulum
double q = 0.5;              // damping coefficient
double Omega_D = 2.0/3.0;    // frequency of driving force
double F_D = 0.9;            // amplitude of driving force
bool nonlinear;              // linear if false

Vector f(const Vector& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    Vector f(3);             // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (nonlinear)
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

Vector x(3);

void getInput() {
    cout << " Nonlinear damped driven pendulum\n"
         << " --------------------------------\n"
         << " Enter linear or nonlinear: ";
    string response;
    cin >> response;
    nonlinear = (response[0] == 'n');
    cout<< " Length of pendulum L: ";
    cin >> L;
    cout<< " Enter damping coefficient q: ";
    cin >> q;
    cout << " Enter driving frequency Omega_D: ";
    cin >> Omega_D;
    cout << " Enter driving amplitude F_D: ";
    cin >> F_D;
    cout << " Enter theta(0) and omega(0): ";
    double theta, omega;
    cin >> theta >> omega;

    x[0] = 0;
    x[1] = theta;
    x[2] = omega;
}

double dt = 0.05;
double accuracy = 1e-6;

void step() {
    adaptiveRK4Step(x, dt, accuracy, f);
}

double frames_per_second = 30;   // for animation in real time

void animation_step() {
    double start = x[0];
    clock_t start_time = clock();
    step();
    double tau = 1.0 / frames_per_second;
    while (x[0] - start < tau)
        step();
    while ((double(clock()) - start_time)/CLOCKS_PER_SEC < tau)
        ;
    glutPostRedisplay();
}

void drawText(const string& str, double x, double y) {
    glRasterPos2d(x, y);
    int len = str.find('\0');
    for (int i = 0; i < len; i++)
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]);
}

void drawVariables() {
    // write t, theta, omega values     
    glColor3ub(0, 0, 255);
    ostringstream os;
    os << "t = " << x[0] << ends;
    drawText(os.str(), -1.1, -1.1);
    os.seekp(0);
    os << "theta = " << x[1] << ends;
    drawText(os.str(), -1.1, 1.1);
    os.seekp(0);
    os << "omega = " << x[2] << ends;
    drawText(os.str(), 0.3, 1.1);
}

void displayPendulum() {
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the pendulum rod
    double xp = sin(x[1]);
    double yp = -cos(x[1]);
    glColor3ub(0, 255, 0);
    glBegin(GL_LINES);
        glVertex2d(0, 0);
        glVertex2d(xp, yp);
    glEnd();

    // draw the pendulum bob
    glPushMatrix();
    glTranslated(sin(x[1]), -cos(x[1]), 0);
    glColor3ub(255, 0, 0);
    const double r = 0.1;
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_POLYGON);
        for (int i = 0; i < 12; i++) {
            double theta = 2 * pi * i / 12.0;
            glVertex2d(r * cos(theta), r * sin(theta));
        }
    glEnd();
    glPopMatrix();

    // write t, theta, and omega
    drawVariables();

    // we are using double buffering - write buffer to screen
    glutSwapBuffers();
}

bool clearPhasePlot;

void displayPhasePlot() {
    static double thetaOld, omegaOld;
    if (clearPhasePlot) {
        glClear(GL_COLOR_BUFFER_BIT);
        clearPhasePlot = false;
        thetaOld = 2 * pi, omegaOld = 2 * pi;
        
        // draw axes
        glColor3ub(0, 255, 0);
        glBegin(GL_LINES);
            glVertex2d(0, -1); glVertex2d(0, 1);
            glVertex2d(-1, 0); glVertex2d(1, 0);
        glEnd();
    }

    // draw the phase point
    double theta = x[1];
    while (theta >= pi) theta -= 2 * pi;
    while (theta < -pi) theta += 2 * pi;
    double omega = x[2];
    glColor3ub(255, 0, 0);
    if (abs(theta - thetaOld) < pi && abs(omega) < pi 
        && abs(omega - omegaOld) < pi) {
        glBegin(GL_LINES);
            glVertex2d(thetaOld / pi, omegaOld / pi);
            glVertex2d(theta / pi, omega / pi);
        glEnd();
    }
    thetaOld = theta, omegaOld = omega;
    glPopMatrix();

    // write t, theta, and omega
    glColor3ub(255, 255, 255);
    glRectd(-1.2, 1, 1.2, 1.2);
    glRectd(-1.2, -1, 1.2, -1.2);
    drawVariables();

    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double aspectRatio = w/double(h);
    double d = 1.2;
    if (aspectRatio > 1)
        glOrtho(-d*aspectRatio, d*aspectRatio, -d, d, -1.0, 1.0);
    else
        glOrtho(-d, d, -d/aspectRatio, d/aspectRatio, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
}

bool running;                    // is the animation running?
bool phasePlot;                  // switch to phase plot if true

void mouse(int button, int state, int x, int y) {
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running) {
                glutIdleFunc(NULL);
                running = false;
            } else {
                glutIdleFunc(animation_step);
                running = true;
            }
        }
        break;
      case GLUT_RIGHT_BUTTON:
        if (state == GLUT_DOWN) {
            if (phasePlot) {
                glutDisplayFunc(displayPendulum);
                phasePlot = false;
            } else {
                glutDisplayFunc(displayPhasePlot);
                clearPhasePlot = phasePlot = true;
            }
            glutPostRedisplay();
        }
        break;
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    getInput();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[0]);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(displayPendulum);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutIdleFunc(NULL);
    glutMainLoop();
}

