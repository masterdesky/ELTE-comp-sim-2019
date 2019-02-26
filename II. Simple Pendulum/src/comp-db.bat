@ECHO OFF
clang-cl -Icpl\\ -Wall pendulum_double.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum_double.exe
@ECHO g++ -Icpl\\ -Wall pendulum_double.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum_double
ECHO Compiling Done!