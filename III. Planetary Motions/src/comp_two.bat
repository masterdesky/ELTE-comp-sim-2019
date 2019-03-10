@ECHO OFF
clang-cl -Icpl\\ -IGL\\ -Wall kepler_two.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler_two.exe
@ECHO g++ -Icpl\\ -Wall kepler_two.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\kepler_two
ECHO Compiling Done!