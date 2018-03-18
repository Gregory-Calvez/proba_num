#pragma once
#include <random>
#include <iostream>
#include <math.h>

#include "random_variable.h"
#include "process.h"
#include "normal.h"

template <typename Generator> class monte_carlo
{
public:
    /// Constructors & destructors
    monte_carlo();
    // monte_carlo(monte_carlo & mc);
    monte_carlo(random_variable<Generator> * rv);
    ~monte_carlo();

    /// Getters & setters
    double get_precision();
    unsigned int get_cap();
    double get_mean();
    double get_std();
    std::pair<double, double> get_confidence_interval();
    void set_precision(double p);
    void set_cap(unsigned int c);

    /// Simulation of one realisation of the random variable
    double simulate(Generator & gen);
    /// Computations of a naive Monte Carlo estimator
    void compute (Generator & gen);
    /// Printing the results
    void print();
public:
    random_variable<Generator> * variable;
    bool computed;
    bool successful;
    double precision;
    unsigned int cap_iterations;
    unsigned int count;
    double empirical_mean;
    double empirical_var;
    double empirical_std;
    double ci_l_bound, ci_h_bound;
};


template <typename Generator> monte_carlo<Generator>::monte_carlo(){
    /// All this are default values
    computed = false;
    successful = false;
    precision = 1.0;
    cap_iterations = 10000;
    count = 0;
    empirical_mean = 0.;
    empirical_var = 0.;
    empirical_std = 0.;
    ci_l_bound = 0.;
    ci_h_bound = 0.;
    /// By default, the variable we intergrate is a N(0, 1)
    variable = new normal<Generator>(0., 1.);
};

template <typename Generator> monte_carlo<Generator>::monte_carlo(random_variable<Generator> * rv){
    /// The random variable we want to integrate
    variable = rv;
    /// All these are default values.
    computed = false;
    successful = false;
    precision = 1.0;
    cap_iterations = 10000;
    empirical_mean = 0.;
    empirical_var = 0.;
    empirical_std = 0.;
    ci_l_bound = 0.;
    ci_h_bound = 0.;
    count = 0;
};

template <typename Generator> monte_carlo<Generator>::~monte_carlo(){};

template <typename Generator> double monte_carlo<Generator>::simulate(Generator & gen){
    return (*variable)(gen);
};

template <typename Generator> void monte_carlo<Generator>::set_precision(double p){
    precision = p;
};

template <typename Generator> void monte_carlo<Generator>::set_cap(unsigned int c){
    cap_iterations = c;
};

template<typename Generator> double monte_carlo<Generator>::get_precision(){
    return precision;
};

template<typename Generator> unsigned int monte_carlo<Generator>::get_cap(){
    return cap_iterations;
};

template<typename Generator> double monte_carlo<Generator>::get_mean(){
    return empirical_mean;
};

template<typename Generator> double monte_carlo<Generator>::get_std(){
    return empirical_std;
};

template<typename Generator> std::pair<double,double> monte_carlo<Generator>::get_confidence_interval(){
    return std::pair<double,double> (ci_l_bound, ci_h_bound);
};

template<typename Generator> void monte_carlo<Generator>::compute(Generator & gen){
    /// This is the method that computes a naive MC estimator
    double sum = 0.;
    double sum_square = 0.;
    double temp = 0.;
    double confidence = 0.;
    do{
        ++count;
        temp = (*variable)(gen);
        sum += temp;
        sum_square += std::pow(temp, 2);
        confidence = 2*1.96*std::sqrt(sum_square/count - std::pow((sum/count), 2)) / std::sqrt(count);
    } while(count < 100 || (confidence > precision && count < cap_iterations));
    empirical_mean = sum / count;
    empirical_var = sum_square / count - std::pow(empirical_mean, 2);
    empirical_std = std::sqrt(empirical_var);
    ci_l_bound = empirical_mean - 1.96 * empirical_std / std::sqrt(count);
    ci_h_bound = empirical_mean + 1.96 * empirical_std / std::sqrt(count);
    successful = (confidence < precision);
    computed = true;
};

template<typename Generator> void monte_carlo<Generator>::print(){
    /// Printing our results properly
    if(!computed){
        std::cout << "Values not computed yet" << std::endl;
    }
    else{
        std::cout << "Values computed, here are the results: " << std::endl;
        if(!successful){
            std::cout << "Precision not reached" << std::endl;
        }
        std::cout << "Empirical mean" << "\t" << this->empirical_mean << std::endl;
        std::cout << "Empirical variance" << "\t" << this->empirical_var << std::endl;
        std::cout << "Number of iterations" << "\t" << this->count << std::endl;
        std::cout << "Confidence interval (95\%)" << "\t (";
        std::cout << this->ci_l_bound << "\t;\t" << this->ci_h_bound << ")" << std::endl;
    }
};
