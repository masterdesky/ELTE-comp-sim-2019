@ECHO OFF
clang-cl -std=c++11 -I..\\..\\cpl\\ -Wall popdin_log.cpp ..\\..\\cpl\\vector.cpp ..\\..\\cpl\\odeint.cpp -o ..\\Release\\popdin_log.exe
clang-cl -std=c++11 -I..\\..\\cpl\\ -Wall popdin_conlog.cpp ..\\..\\cpl\\vector.cpp ..\\..\\cpl\\odeint.cpp -o ..\\Release\\popdin_conlog.exe
clang-cl -std=c++11 -I..\\..\\cpl\\ -Wall popdin_lv.cpp ..\\..\\cpl\\vector.cpp ..\\..\\cpl\\odeint.cpp -o ..\\Release\\popdin_lv.exe
ECHO Compiling Done!