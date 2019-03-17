#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>

#include "vector.hpp"
#include "odeint.hpp"

std::vector<cpl::Vector> SmallBodies(100, cpl::Vector(5));
std::vector<double> SmallBodyMasses(100);

int main(){
    // Read in data of small bodies into a vector
    std::ifstream inputFile("small_objects.dat");

    int number_of_bodies = 0;                   // Number of small bodies in the simulation
    // Check if input file exists
    if(inputFile.good()) {
        
        double current_number = 0;              // Storage for current read number from data file
        int current_column = 0;
        int current_row = 0;
        // Push items into a vector
        while (inputFile >> current_number) {
            std::cout << current_number << '\t';

            if(current_column == 0) {
                SmallBodyMasses[current_row] = current_number;
                SmallBodies[current_row][0] = 0;
            }

            else{
                SmallBodies[current_row][current_column] = current_number;
            }

            current_column += 1;

            if(current_column == 5) {
                current_column = 0;
                current_row++;
                number_of_bodies++;
                std::cout << '\n';
            }
        }

        // Close the file.
        inputFile.close();
    }

    else {
        exit(-1);
    }

    // Print read in datapoints
    std::cout << "Read in data:\n"
              << " -------------------------------\n";
    for(int i = 0; i < number_of_bodies; i++) {
        for(int j = 0; j < 5; j++) {
            std::cout << SmallBodies[i][j] << '\t';
        }
        std::cout << '\n';
    }
}