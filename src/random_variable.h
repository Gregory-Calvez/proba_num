#pragma once
#include <random>
#include <iostream>
#include <math.h>
#include <utility>

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
    virtual std::pair<double, double> simulate(Generator & gen);

    double get_mean_control_variate();
    bool has_control();
protected:
    bool realised;
    double realisation;
    bool has_control_variate;

    // For control variate
    double control_variate_realisation;
    double control_variate_mean;
};




template <typename Generator> random_variable<Generator>::~random_variable() {

};

template <typename Generator> random_variable<Generator>::random_variable() {
    this->realised = false;
    this->realisation = 0;
    this->has_control_variate = false;
    this->control_variate_mean = 0;
    this->control_variate_realisation = 0;
};

template <typename Generator> void random_variable<Generator>::print(){
    if(realised) {
        std::cout << realisation << std::endl;
    };
};

template <typename Generator> bool random_variable<Generator>::has_control(){
    return this->has_control_variate;
};

template <typename Generator> double random_variable<Generator>::get_mean_control_variate(){
    return this->control_variate_mean;
};

template <typename Generator> std::pair<double, double> random_variable<Generator>::simulate(Generator & gen){
    if(this->has_control()){
        std::cout << "This should be overload" << std::endl;
        return std::pair<double, double> ((*this)(gen), 0);
    }
    else {
        return std::pair<double, double> ( (*this)(gen), 0);
    };
}
