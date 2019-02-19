@ECHO OFF
@ECHO clang-cl -Icpl\ -Wall -o ..\Release\pendulum.exe pendulum.cpp
g++ -Icpl\\ -Wall pendulum.cpp -o ..\\Release\\pendulum
ECHO Compiling Done!