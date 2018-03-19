#pragma once
#include <random>
#include <iostream>
#include <math.h>

#include "random_variable.h"
#include "monte_carlo.h"
#include "process.h"
#include "brownian.h"

/** \brief A class to simulate the integral of a Brownian Motion with a classical Euler scheme.
*
* Mainly used for testing.
*/
template<typename Generator> class integral_brownian: public random_variable<Generator>
{
/*
This class generates simulation of the integral b/w 0 and 1 of a brownian motion.
It uses the class brownian.
This is a class inherited from random_variable.
It is used for monte_carlo in the main function.
*/
public:
    integral_brownian();
    integral_brownian(unsigned int number_steps);
    ~integral_brownian();
    // integral_brownian(integral_brownian<Generator> & ib;

    double operator() (Generator & gen);
private:
    /// With default t_max = 1
    brownian<Generator> * b;
    unsigned int num_steps;
    double time_step;
};



template<typename Generator> integral_brownian<Generator>::integral_brownian(unsigned int number_steps){
    this->num_steps = number_steps;
    this->b = new brownian<Generator>();
    b->set_num_steps(this->num_steps);
    time_step = b->get_time_step();
};

template<typename Generator> integral_brownian<Generator>::integral_brownian(){
    this->integral_brownian(100);
};

template<typename Generator> integral_brownian<Generator>::~integral_brownian(){
    delete(b);
};

template<typename Generator> double integral_brownian<Generator>::operator() (Generator & gen){
    std::vector<std::vector<double> > x = (*(this->b))(gen);
    double s = 0;
    for(std::vector<std::vector<double> >::iterator it = x.begin(); it != x.end(); ++it){
        s += (*it).at(0)*this->time_step;
    };
    this->realisation = s;
    this->realised = true;
    return s;
};
