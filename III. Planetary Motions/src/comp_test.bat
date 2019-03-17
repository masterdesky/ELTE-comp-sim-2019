@ECHO OFF
clang-cl -std=c++11 -Icpl\\ -IGL\\ -Wall read_test.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\read_test.exe
@ECHO g++ -Icpl\\ -Wall read_test.cpp cpl\\vector.cpp cpl\\odeint.cpp -o ..\\Release\\read_test
ECHO Compiling Done!