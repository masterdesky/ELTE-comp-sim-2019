@ECHO OFF
clang-cl -std=c++11 -Wall cell.cpp -o ..\\Release\\cell.exe
clang-cl -std=c++11 -Wall sandpile.cpp -o ..\\Release\\sandpile.exe
ECHO Compiling Done!