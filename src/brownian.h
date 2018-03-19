#pragma once
#include "random_variable.h"
#include "process.h"

/** \brief A class to discretize a Brownian Motion with a classical Euler scheme.
*
* Mainly used for testing. 
*/

template<typename Generator> class brownian : public process<Generator>
{
public:
    /// Constructors & destructors
    brownian();
    ~brownian();
    // brownian(brownian<Generator> & b);
    brownian(double vol);

    /// Next step
    std::vector<double> next_step(Generator & gen, std::vector<double> state, double t);
private:
    std::normal_distribution<double> increment;
    double volatility;
};

////////////////////////////////////
/// Functions for the brownian class
///////////////////////////////////
template<typename Generator> brownian<Generator>::brownian() : brownian<Generator> (1.){

};
template<typename Generator> brownian<Generator>::~brownian(){
};
template<typename Generator> brownian<Generator>::brownian(double vol){
    volatility = vol;
    increment = std::normal_distribution<double>(0.0, 1.0);
};
const std::vector<double> vector_plus(std::vector<double> v_1, std::vector<double> v_2){
    unsigned int n_1 = v_1.size();
    unsigned int n_2 = v_2.size();
    if(n_1 != n_2){
        return std::vector<double>(0.);
    }
    else {
        std::vector<double> result = std::vector<double> ();
        for(unsigned int i = 0; i<n_1; ++i){
            result.push_back(v_1.at(i) + v_2.at(i));
        };
        return result;
    };
}
template<typename Generator> std::vector<double> brownian<Generator>::next_step(Generator & gen, std::vector<double> state, double t){
    std::vector<double> delta = std::vector<double> (this->dimension, this->volatility * std::sqrt(this->time_step) * increment(gen));
    std::vector<double> result = vector_plus(state, delta);
    // std::cout << result.at(0) << std::endl;
    return result;
};
