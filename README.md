# Computer Simulations 2018-2019/2 @ ELTE
#### Language of course and documentations: Hungarian

The course's initial purpose is to get its students acquainted with numerical simulations of classical problems in physics. During the course, one needs to solve various tasks from 6 different topics, based on a simulation source code for every topic, which was created and made available in advance by the course's tutors. One needs to alter the source code in numerous ways, and analyze the output data to solve the given objectives. It is also required to make a documentation about the whole process for every subjects.  
The titles and abstracts of the topics are listed below. **For the original `C++` source codes, all credit goes to the tutors of the course!**

## I. Simulation: Harmonic Oscillator
On the first class of the *Computer Simulations* course we studied the simple 1D harmonic oscillator and the numerical solution of its equation of motion. The class serves as an introductory for the course. The source code for the simulation was written on `C++` language, and the data analysis and the plotting was made with a Python 3 kernel running inside Jupyter Notebook. The first submission initial task was to get well known with the source code and then excercise by changing its input variables to achieve various runs with different initial conditions. Following this, we was needed to analyse the outputs according to given aspects and conditions.

## II. Simulation: Simple Pendulum
On the second class of the *Computer Simulations* course we examined the numerical solutions of the differential equations of simple pendulums. We studied various approximations of the subject, like simple gravity pendulum, driven and/or damped pendulum and physical pendulum. We examined in contrast the simple fourth-order Runge-Kutta and stepsize controlled (adaptive) fourth-order Runge-Kutta, the Runge-Kutta-Cash-Karp and adaptive Runge-Kutta-Cash-Karp and the Euler and Euler-Cromer algorithms. Furthermore we also reviewed the same algorithms for the double pendulum.

## III. Simulation: Planetary Motions
On the third occasion of the *Computer Simulations* course, numerical methods for simple planetary motion and the Kepler problem were investigated. We have compared the stepsize controlled method - that has been used in the previous simulations - to solve the differential equations of these problems and test, how it behaves under different initial parameters. In addition, we tested the runtime in contrast with the fixed step integration method.  
We modeled the one-, two-, and three-body problems, as well as examined Mercury's perihelion precession, supplementing the motion equations with relativistic effects. Finally, in the case of the three-body problem, we modeled the movement of the asteroids and meteorites in the Sun-Jupiter system and monitored their equilibrium at the system's Lagrange points.

## IV. Simulation: Molecular Dynamics
On the fourth occasion of the *Computer Simulations* laboratory we reviewed thermodynamics and statistical physics, and their approach to physics. In this case, this meant that we studied the microscopic dynamics of large number of interacting particles and used it to map some macroscopic properties of the examined system. These included examining the equilibrium position of the system and modeling temperature, pressure, total energy, compression factor, and heat capacity in this state. In addition, we got to know the Verlet and velocity-Verlet algorithms, as well as the corrections used along with them. These algorithms are commonly used in molecular dynamics simulations because of their speed, compared to other iterative methods.

## V. Simulation: Population Dynamics
On the fifth occasion of the *Computer Simulations* course we dealt with a topic outside physics, which is part of the science of biology. This is the topic of population dynamics, which deals with the temporal changes of individual animal (or human) populations, and their interaction with the environment and other species. The topic has gained a reason to get into the subject itself by the fact, that such temporal changes can easily be modeled with numerical solutions of differential equation. During the tasks, we also examined the logistic, the attached-logistic and the Lotka-Volterra model, which are the most basic population dynamics algorithms.

## VI. Simulation: Cellular Automaton
On the sixth and last occasion of the *Computer Simulations* course we examined simple, 2D cellular automatons. We've implemented Conway's Game of Life and a sand pile automaton.
