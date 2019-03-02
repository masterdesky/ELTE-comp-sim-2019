@ECHO OFF
clang-cl -Icpl\\ -IGL\\ -Wall kepler.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler.exe
@ECHO g++ -Icpl\\ -Wall kepler.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler
ECHO Compiling Done!