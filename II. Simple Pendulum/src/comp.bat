@ECHO OFF
clang-cl -Icpl\\ -IGL\\ -Wall pendulum.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum.exe
@ECHO g++ -Icpl\\ -Wall pendulum.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum
ECHO Compiling Done!