@ECHO OFF
clang-cl -std=c++11 -Icpl\\ -IGL\\ -Wall kepler_three.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler_three.exe
@ECHO g++ -Icpl\\ -Wall kepler_three.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler_three
ECHO Compiling Done!