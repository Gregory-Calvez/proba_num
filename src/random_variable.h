#pragma once
#include <random>
#include <iostream>
#include <math.h>

/** \brief A general class to implement a random variable
*/

template<typename Generator> class random_variable
{
public:
    /// Constructor
    random_variable();
    /// Destructor
    ~random_variable();
    /// Printing
    void print();
    /// Virtual function
    virtual double operator()(Generator & gen) = 0;
protected:
    bool realised;
    double realisation;
};




template <typename Generator> random_variable<Generator>::~random_variable() {

};

template <typename Generator> random_variable<Generator>::random_variable() {
    this->realised = false;
    this->realisation = 0;
};

template <typename Generator> void random_variable<Generator>::print(){
    if(realised) {
        std::cout << realisation << std::endl;
    };
};
