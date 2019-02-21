#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
using namespace std;

#include "random.hpp"

namespace cpl {

int timeSeed() {
    time_t t;                   // usually an unsigned long int
    time(&t);                   // get seconds since Jan 1, 1900
    tm* pt = gmtime(&t);        // Universal (Greenwich Mean) Time
    return pt->tm_year + 70 * (pt->tm_mon + 12 * (pt->tm_hour +
           23 * (pt->tm_min + 59 * pt->tm_sec)));
}

double rand(int& seed, bool set) {
    if (set)
        srand(seed);
    else
        seed = std::rand();
    return seed / (RAND_MAX + 1.0);
}

double ranf(int& seed, bool set) {
    const int IA = 16807, IC = 2147483647, IQ = 127773, IR = 2836;
    if (!set) {
        int h = seed / IQ;
        int l = seed % IQ;
        seed = IA * l - IR * h;
    }
    if (seed <= 0)
        seed += IC;
    return seed / double(IC);
}

double qand(int& seed, bool set) {
    static unsigned long idum;
    const double TWO_POWER_32 = 4294967296.0;
    if (set) {
        idum = (unsigned long) seed;
        return idum / TWO_POWER_32;
    }
    idum = 1664525L * idum + 1013904223L;
    seed = int(idum);
    return idum / TWO_POWER_32;
}

double gasdev(int& seed) {
     static int iset = 0;
     static double gset;
     double fac, rsq, v1, v2;
     if (iset == 0) {
          do {
               v1 = 2.0 * ranf(seed) - 1.0;
               v2 = 2.0 * ranf(seed) - 1.0;
               rsq = v1 * v1 + v2 * v2;
          } while (rsq >= 1.0 || rsq == 0.0);
          fac = sqrt(-2.0 * log(rsq)/rsq);
          gset = v1 * fac;
          iset = 1;
          return v2 * fac;
     } else {
          iset = 0;
          return gset;
     }
}

}  /* end namespace cpl */
