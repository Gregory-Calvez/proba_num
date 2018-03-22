#pragma once
#include <random>
#include <iostream>
#include <math.h>

#include "random_variable.h"
#include "process.h"
#include "normal.h"

/** \brief A general class for Monte-Carlo simulation
*
* In this class we compute naive MC estimators (with confidence interval) given a random variable.
* We also implement different variance reduction techniques.
*/

template <typename Generator> class monte_carlo
{
public:
    /// Constructors & destructors
    monte_carlo();
    // monte_carlo(monte_carlo & mc);
    monte_carlo(random_variable<Generator> * rv);
    ~monte_carlo();

    /// Getter
    double get_precision();
    /// Getter
    unsigned int get_cap();
    /// Getter
    double get_mean();
    /// Getter
    double get_std();
    /// Getter
    std::pair<double, double> get_confidence_interval();
    /// Setter
    void set_precision(double p);
    /// Setter
    void set_cap(unsigned int c);

    /// Simulation of one realisation of the random variable
    double simulate(Generator & gen);
    /// Computations of a naive Monte Carlo estimator
    void compute (Generator & gen);
    /// Computations of a Monte Carlo with control variate
    void compute_control_variate (Generator & gen);
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
    /// By default, the variable we integrate is a N(0, 1)
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


// To simplify things, we will say that a random variable can now return a pair with
// as first value the reaslisation of the RV and as second value the value of
// the control variate.
// We will say that we have only one control variate.
template<typename Generator> void monte_carlo<Generator>::compute_control_variate(Generator & gen)
{
    if (this->variable->has_control() == false){
        std::cout << "This variable has no control variate defined. " << std::endl;
        return;
    };

    // First, compute the control constant
    double mean_control_variate = this->variable->get_mean_control_variate();
    unsigned int num_iterations_to_compute_constant = std::floor(this->cap_iterations / 10);

    double sum_variable = 0;
    double sum_control = 0;
    double sum_covariance = 0;
    double sum_order_2_control = 0;
    std::pair<double, double> temp;
    for(unsigned int i = 0; i<num_iterations_to_compute_constant; ++i){
        temp = this->variable->simulate(gen);
        sum_variable += temp.first;
        sum_control += temp.second;
        sum_covariance += temp.first*temp.second;
        sum_order_2_control += std::pow(temp.second, 2);
    };
    double covariance = (sum_covariance / num_iterations_to_compute_constant) - mean_control_variate * (sum_variable / num_iterations_to_compute_constant);
    double variance_control = (sum_order_2_control / num_iterations_to_compute_constant) - std::pow(mean_control_variate, 2);

    double c = - covariance / variance_control;


    // Then, we do the monte carlo
    double sum_mc = 0.;
    double sum_mc_square = 0.;
    double confidence = 0.;
    double temp_x = 0;
    unsigned int count = 0;
    do{
        ++count;
        temp = this->variable->simulate(gen);
        temp_x = temp.first + c * (temp.second - mean_control_variate);
        sum_mc += temp_x;
        sum_mc_square += std::pow(temp_x, 2);
        confidence = 2*1.96*std::sqrt(sum_mc_square/count - std::pow((sum_mc/count), 2)) / std::sqrt(count);
    } while(count < 100 || (confidence > precision && count < cap_iterations));
    empirical_mean = sum_mc / count;
    empirical_var = sum_mc_square / count - std::pow(empirical_mean, 2);
    empirical_std = std::sqrt(empirical_var);
    ci_l_bound = empirical_mean - 1.96 * empirical_std / std::sqrt(count);
    ci_h_bound = empirical_mean + 1.96 * empirical_std / std::sqrt(count);
    successful = (confidence < precision);
    computed = true;
    return ;
};
