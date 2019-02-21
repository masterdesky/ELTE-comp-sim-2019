#ifndef CPL_VECTOR_HPP
#define CPL_VECTOR_HPP

#include <iostream>

namespace cpl {

// Vector object with components of type double

class Vector {
  public:

    Vector(                // constructor, creates a new vector
        int dim = 1        // dimension, assumed = 1 if omitted
    );                     // components initialized to 0

    Vector(                // constructor, creates a new 2-d vector
        double v0,         // value of first component
        double v1          // value of second componet
    );

    Vector(                // constructor, creates a new 3-d vector
        double v0,         // value of first component
        double v1,         // value of second component
        double v2          // value of third component
    );

    Vector(                // copy constructor, creates a new vector
        const Vector&      // vector to be copied
    );

    ~Vector() {            // destructor frees memory allocated by constructors
        delete [] v;       // inline code to delete array v
    }

    int dimension() const  // returns number of components
    { return dim; }        // inline code to get dimension

    int size() const       // same as dimension()
    { return dim; }

    void resize(           // resize the vector
        const int          // new number of components
    );                     // old components preserved, any new zeroed

    void push_back(        // add a new component to the vector
        const double       // value to set the new component
    );

    void set(              // set all the components of the vector
        const double,      // value of first component
        ...                // second, third, ... , last component
    );                     // all components must be specified

    const double           // value returned by
    operator[](            // [ ] operator function
        const int i        // index of component to return
    ) const                // value cannot be modified
    { return v[i]; }       // inline code to get value

    double&                // reference returned by
    operator[](            // [ ] operator function
        const int i        // index of component to return
    ) { return v[i]; }     // inline code to get component

    // operators to allow the vector to be assigned

    Vector& operator = (const Vector&);

    Vector& operator += (const Vector&);

    Vector& operator -= (const Vector&);

    // multiply or divide the vector by a scalar

    Vector& operator *= (double);

    Vector& operator /= (double);

    double abs();          // return the length of the vector

    double norm();         // return the norm of the vector

    double dot(            // return dot product of the vector
       const Vector&       // with this vector
    );

    // allow the << operator to be used to output vector components
    friend std::ostream& operator<<(std::ostream& os, const Vector& dv);

  private:
    int dim;               // store the number of components of the vector
    double *v;             // pointer to stored components values
};

// various arithmetic operations on vectors

// unary plus using inline code, just returns the same vector
inline Vector operator + (const Vector& vec) { return vec; }

// unary minus, flips the sign of each component
extern Vector operator - (const Vector&);

// multiplying and dividing a vector by a scalar

extern Vector operator * (const Vector&v, double);

extern Vector operator * (double, const Vector&);

extern Vector operator / (const Vector&, double);

// adding and subtracting two vectors

extern Vector operator + (const Vector&, const Vector&);

extern Vector operator - (const Vector&, const Vector&);

// dot product of two vectors

extern double dot(const Vector&, const Vector&);

} /* end namespace cpl */

#endif /* CPL_VECTOR_HPP */
