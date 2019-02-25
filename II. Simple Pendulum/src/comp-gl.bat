@ECHO OFF
clang-cl -Icpl\\ -IGL\\ -Wall pendulum-gl.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum-gl.exe
@ECHO g++ -Icpl\\ -IGL\\ -Wall pendulum-gl.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\pendulum-gl
ECHO Compiling Done!