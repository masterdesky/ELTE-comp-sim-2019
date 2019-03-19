@ECHO OFF
clang-cl -std=c++11 -Icpl\\ -IGL\\ -Wall kepler.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler.exe
ECHO Compiling Done!