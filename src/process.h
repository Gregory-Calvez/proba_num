#pragma once
#include <random>
#include <iostream>
#include <math.h>

#include "random_variable.h"
#include "monte_carlo.h"

/** \brief A general class to code the discretization of a stochastic process on the time interval [0,T]
*
* The discretization scheme can be fixed either by chosing the time horizon (time_end)
* or the number of steps (num_steps).
*/


template<typename Generator>
class process
{
public:
    /// Constructor
    process();
    ///Destructor
    ~process();


    /// Setter
    void set_time_end(double t);
    /// Setter
    void set_num_steps(unsigned int n);
    /// Setter
    void set_state_0(std::vector<double> s);
    /// Setter
    void set_time_grid();
    /// Setter
    void set_state(std::vector<double> s);
    /// Getter
    unsigned int get_num_steps();
    /// Getter
    std::vector<double> get_state_0();
    /// Getter
    double get_time_end();
    /// Getter
    double get_time_step();

    //! For plotting
    /*!
      Computes and plots num_plot instances of the discretized process.
      The method writes the relevant data in a .dat file and calls gnuplot for plotting
    */
    void write_to_plot(Generator & gen, unsigned int num_plots, int coordinate);
    //! Function Next Step
    /*!
    *  Pure virtual function that will be over-written for every instance of process.
    *  Computes the value of the discretized process at time t + time_step.
    *  The value of the discretized process at t is stored in state.
    */
    virtual std::vector<double> next_step(Generator & gen, std::vector<double> state, double t) = 0;
    //! Operator()
    /*!
      Overloads the operator ().
      Allows to draw a new instance of the discretized process.
    */
    std::vector<std::vector<double> > operator() (Generator & gen);
    //! Prints the results of the discretization
    void print_trajectory(int coordinate);

protected:
    unsigned int dimension;
    double time_step; /**< Time step of the discretization */
    double time_end; /**<Simulation between 0 and time_end */
    unsigned int num_steps; /**<Number of steps in time discretization*/
    std::vector<double> state_0; /**<Initial value for the process at time t = 0 */
    std::vector<double> state;
    std::vector<double> time_grid; // time grid
    std::vector<std::vector<double> > trajectory; /**<Vector of vectors used to compute the values taken by the process
                                                  This choice is motivated by the mutli-dimensional Heston model */
    bool computed; // boolean: true if trajectory has been computed
};

template<typename Generator> void process<Generator>::set_state(std::vector<double> s)
{
    this->state = s;
};
template<typename Generator> void process<Generator>::write_to_plot(Generator & gen, unsigned int num_plots, int coordinate)
{
    std::vector<std::vector<std::vector<double> > > data = std::vector<std::vector<std::vector<double > > > ();
    for(unsigned int i = 0; i<num_plots; ++i){
        (*this)(gen);
        data.push_back(trajectory);
    };
    std::cout << "We have computed the trajectories. " << std::endl;

    std::ofstream stream;
    stream.open("plot.dat");
    stream << "time";
    for(unsigned int i = 0; i<num_plots; ++i){
        stream << "\t" << (i+1);
    };
    stream << "\n";
    for (unsigned int i_t = 0; i_t < this->num_steps; ++i_t){
        stream << this->time_grid.at(i_t);
        for(unsigned int i = 0; i<num_plots; ++i){
            stream << "\t" << data.at(i).at(i_t).at(coordinate);
        };
        stream << "\n";
    };
    stream.close();
    std::cout << "We have written the data. " <<std::endl;
    stream.open("gnu");
    stream << "set nokey\n" ;
    stream << "set xlabel \"Time\"\n";
    stream << "plot ";
    for(unsigned int i = 1; i<num_plots+1; ++i){
        stream << "\"plot.dat\" using 1:" << i+1;
        stream << " with lines, \\\n";
    };
    stream.close();
    std::cout << "We have written the gnuplot file. " << std::endl;
};



template<typename Generator> void process<Generator>::set_time_grid(){
    time_grid = std::vector<double>();
    trajectory = std::vector<std::vector<double> > ();
    double time = 0.;
    for(unsigned int i = 0; i<this->num_steps; ++i){
        time_grid.push_back(time);
        time += time_step;
        trajectory.push_back(state_0);
    };
};
template<typename Generator> void process<Generator>::set_time_end(double t){
    time_end = t;
    time_step = time_end / num_steps;
    this->set_time_grid();
};
template<typename Generator> void process<Generator>::set_num_steps(unsigned int n){
    this->num_steps = n;
    time_step = time_end / num_steps;
    this->set_time_grid();
};
template<typename Generator> void process<Generator>::set_state_0(std::vector<double> s){
    state_0 = s;
};

template<typename Generator> unsigned int process<Generator>::get_num_steps(){
    return num_steps;
};

template<typename Generator> double process<Generator>::get_time_end(){
    return time_end;
};

template<typename Generator> double process<Generator>::get_time_step(){
    return time_step;
};

template<typename Generator> std::vector<double> process<Generator>::get_state_0(){
    return state_0;
};

template<typename Generator> void process<Generator>::print_trajectory(int coordinate){
    if(computed){
        std::cout << "Values of the process" << std::endl;
        std::cout << "t\tprocess" << std::endl;
        for(unsigned int i = 0; i < num_steps; ++i){
            std::cout << time_grid.at(i) << "\t" << trajectory.at(i).at(coordinate) << std::endl;
        };
    }
    else{
        std::cout << "Values not computed yet" << std::endl;
    };
};
template<typename Generator> process<Generator>::process(){
    this->dimension = 1;
    time_end = 1.0;
    num_steps = 100;
    time_step = .01;
    time_grid = std::vector<double>();
    trajectory = std::vector<std::vector<double> >();
    state_0 = std::vector<double> (1.);
    double time = 0.;
    for(unsigned int i = 0; i<num_steps; ++i){
        time_grid.push_back(time);
        time += time_step;
        trajectory.push_back(state_0);
    };
    computed = false;
    this->state = this->state_0;
};
template<typename Generator> process<Generator>::~process(){
};
template<typename Generator> std::vector<std::vector<double> > process<Generator>::operator() (Generator & gen){
    trajectory = std::vector<std::vector<double> >();
    this->state = state_0;
    double t = 0;
    for(unsigned int i = 0; i<num_steps; ++i){
        trajectory.push_back(this->state);
        this->state = this->next_step(gen, this->state, t);
        t += time_step;
    };
    trajectory.push_back(this->state);
    computed = true;
    return trajectory;
};
