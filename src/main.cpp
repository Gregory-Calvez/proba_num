#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <random>


#include "normal.h"
#include "random_variable.h"
#include "monte_carlo.h"
#include "process.h"
#include "brownian.h"
#include "cir.h"
#include "heston.h"
#include "integral_brownian.h"
#include "option.h"
#include "sobol_generator.h"

void plot_heston_o2(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, unsigned int num_steps, int coordinate){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Construct Heston
    heston<std::mt19937_64, cir_o2<std::mt19937_64> > h = heston<std::mt19937_64, cir_o2<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    h.set_num_steps(num_steps);
    h(generator);
    // Generate and write the plots
    // 0 = CIR, 1 = \int CIR, 2 = Heston, 3 = \int Heston
    h.write_to_plot(generator, 10, coordinate);
};

void plot_heston_o3(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, unsigned int num_steps, int coordinate){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Construct Heston
    heston<std::mt19937_64, cir_o3<std::mt19937_64> > h = heston<std::mt19937_64, cir_o3<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    h.set_num_steps(num_steps);
    h(generator);
    // Generate and write the plots
    // 0 = CIR, 1 = \int CIR, 2 = Heston, 3 = \int Heston
    h.write_to_plot(generator, 10, coordinate);
};

void plot_heston_glasserman(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, unsigned int num_steps, int coordinate){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Construct Heston
    heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > h = heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    h.set_num_steps(num_steps);
    h(generator);
    // Generate and write the plots
    // 0 = CIR, 1 = \int CIR, 2 = Heston, 3 = \int Heston
    h.write_to_plot(generator, 10, coordinate);
};

void compare_reduction_variance_glasserman(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, double strike, double expiry, unsigned int num_steps, char type, unsigned int cap){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Initialize the option
    option<std::mt19937_64, heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > > * opt = new option<std::mt19937_64, heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    monte_carlo<std::mt19937_64> mc(opt);
    mc.set_precision(0.01);
    mc.set_cap(cap);
    // Run without control variates
    std::cout << "Without control variates" << std::endl;
    mc.compute(generator);
    mc.print();
    std::cout << std::endl << "With control variates";
    // Run with control variates

    mc.compute_control_variate(generator);
    mc.print();
    std::cout << std::endl;
};

void compare_reduction_variance_o2(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, double strike, double expiry, unsigned int num_steps, char type, unsigned int cap){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Initialize the option
    option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt = new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    if (type == 'e'){
        std::cout << "Test of Monte-Carlo on an European option" << std::endl;
    } else {
        std::cout << "Test of Monte-Carlo on an Asian option" << std::endl;
    }
    monte_carlo<std::mt19937_64> mc(opt);
    mc.set_precision(0.01);
    mc.set_cap(cap);
    // Run without control variates
    std::cout << "Without control variates" << std::endl;
    mc.compute(generator);
    mc.print();
    std::cout << std::endl << "With control variates";
    // Run with control variates

    mc.compute_control_variate(generator);
    mc.print();
    std::cout << std::endl;
};

void compare_reduction_variance_o3(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, double strike, double expiry, unsigned int num_steps, char type, unsigned int cap){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Initialize the option
    option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > * opt = new option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    if (type == 'e'){
        std::cout << "Test of Monte-Carlo on an European option" << std::endl;
    } else {
        std::cout << "Test of Monte-Carlo on an Asian option" << std::endl;
    }
    monte_carlo<std::mt19937_64> mc(opt);
    mc.set_precision(0.01);
    mc.set_cap(cap);
    // Run without control variates
    std::cout << "Without control variates" << std::endl;
    mc.compute(generator);
    mc.print();
    std::cout << std::endl << "With control variates" << std::endl;
    // Run with control variates
    mc.compute_control_variate(generator);
    mc.print();
    std::cout << std::endl;
};

void test_for_cuda(){
    double k = 0.5;
    double a = 0.02;
    double sigma = 0.4;
    double x_0 = 100.;
    double cir_0 = .04;
    double rho = -0.5;
    double r = 0.02;
    double strike = 100.;
    double expiry = 1.;
    unsigned int cap = 2621440;

    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);
    // Initialize the option
    option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt = new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, 'e');
    monte_carlo<std::mt19937_64> mc(opt);
    mc.set_precision(0.0001); // To make sure that we hit the cap for the number of iteration to compare with the other implementations.
    mc.set_cap(cap);
    mc.compute(generator);
    mc.print();
};

void test_mc_sobol(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, double strike, double expiry, unsigned int num_steps, unsigned int cap){
    // Initialize generator, random seed
    sobol_generator sobol = sobol_generator();
    // Construct Heston
    option<sobol_generator, heston<sobol_generator, cir_o2<sobol_generator> > > * opt = new option<sobol_generator, heston<sobol_generator, cir_o2<sobol_generator> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, 'e');
    heston<std::mt19937_64, cir_o2<std::mt19937_64> > h = heston<std::mt19937_64, cir_o2<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    h.set_num_steps(num_steps);
    opt->set_num_steps(num_steps);
    monte_carlo<sobol_generator> mc(opt);
    mc.set_precision(0.0001);
    mc.set_cap(cap);
    std::cout << "Sobol Monte Carlo on CIR_O2 Heston Asian" << std::endl;
    mc.compute(sobol);
    mc.print();
};

void plot_graph_performance(double expiry, double strike, double cir_0, double x_0, double a, double k , double sigma, double rho,double r, char type, std::vector<int> num_steps, unsigned int cap, double precision, double exact_value){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);

    std::vector<std::vector<double> > data = std::vector<std::vector<double> > ();

    option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > * opt_o3 =
        new option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type);

    option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt_o2 =
        new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type);


    for (unsigned int i = 0; i < num_steps.size(); ++i){
        //CIR_O2
        opt_o2->set_num_steps(num_steps.at(i));
        monte_carlo<std::mt19937_64> mc(opt_o2);
        mc.set_precision(precision);
        mc.set_cap(cap);
        mc.compute_control_variate(generator);
        std::vector<double> temp = std::vector<double>();
        temp.push_back(mc.empirical_mean);
        temp.push_back(mc.empirical_mean - mc.ci_l_bound);

        //CIR_03
        opt_o3->set_num_steps(num_steps.at(i));
        monte_carlo<std::mt19937_64> mc_3(opt_o3);
        mc_3.set_precision(precision);
        mc_3.set_cap(cap);

        mc_3.compute_control_variate(generator);

        temp.push_back(mc_3.empirical_mean);
        temp.push_back(mc_3.empirical_mean - mc_3.ci_l_bound);

        data.push_back(temp);
    }
    std::cout << "We have computed the trajectories. " << std::endl;

    /// Writing the data into the plot file.
    std::ofstream stream;
    stream.open("plot.dat");
    for (unsigned int i_n = 0; i_n < num_steps.size(); ++i_n ){
        stream << 1./num_steps.at(i_n);
        stream << "\t" << data.at(i_n).at(0);
        stream << "\t" << data.at(i_n).at(1);
        stream << "\t" << data.at(i_n).at(2);
        stream << "\t" << data.at(i_n).at(3);
        stream << "\n";
    }
    stream.close();
    std::cout << "We have written the data. " <<std::endl;
    stream.open("gnu");
    stream << "set nokey\n" ;
    stream << "set xlabel \"Inverse of number of steps\"\n";
    stream << "plot ";
    stream << "\"plot.dat\" using 1:2:3 with yerrorlines, \\\n";
    stream << "\"plot.dat\" using 1:4:5 with yerrorlines, \\\n";
    stream << exact_value << "with lines lt 3";
    stream.close();
    std::cout << "We have written the gnuplot file. " << std::endl;
};


void pricing_asian(double cir_0, double x_0, double a, double k, double sigma, double rho, double r, double strike, double expiry, unsigned int num_steps, unsigned int cap){
    // Initialize generator, random seed
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);

    // Construct Heston
    option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt = new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (expiry, strike, cir_0, x_0, a, k , sigma, rho, r, 'a');
    heston<std::mt19937_64, cir_o2<std::mt19937_64> > h = heston<std::mt19937_64, cir_o2<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    h.set_num_steps(num_steps);
    opt->set_num_steps(num_steps);
    monte_carlo<std::mt19937_64> mc(opt);
    mc.set_precision(0.0001);
    mc.set_cap(cap);
    std::cout << "Pricing an Asian option. " << std::endl;
    mc.compute(generator);
    mc.print();
};

int main(){
    /// Plotting graphs to measure the performances of the different schemes
    /// Let's define the parameters
    double k = 0.5;
    double a = 0.02;
    double sigma = 0.4;
    double x_0 = 100.;
    double cir_0 = .04;
    double rho = -0.5;
    double r = 0.02;
    double strike = 100.;
    double expiry = 1.;
    unsigned int num_steps = 100;
    unsigned int cap = 100000;


    /// Plotting some examples of trajectories.
    /// The last parameter is the coordinate : 0 = CIR, 2 = Heston
    // std::cout << "This plots some examples of trajectories for the Heston/Cir processes" << std::endl;
    // plot_heston_glasserman(cir_0, x_0, a, k, sigma, rho, r, num_steps, 2);
    // plot_heston_o2(cir_0, x_0, a, k, sigma, rho, r, num_steps, 0);
    // plot_heston_o3(cir_0, x_0, a, k, sigma, rho, r, num_steps, 2);


    /// Testing the control variates
    // std::cout << "This compares the variances for with/out control variates for two options." << std::endl;
    // compare_reduction_variance_glasserman(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'e', cap);
    // compare_reduction_variance_glasserman(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'a', cap);

    // compare_reduction_variance_o2(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'e', cap);
    // compare_reduction_variance_o3(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'a', cap);

    /// One shot to compare with cuda code
    // test_for_cuda();

    /// Test for Sobol
    // std::cout << "This runs one monte carlo with the sobol sequence as RNG" << std::endl;
    // cap  = 10000;
    // test_mc_sobol(cir_0, x_0, a, k, sigma, rho, r, strike,  expiry, num_steps, cap);


    /// Test for Alfonsi's graphs
    // std::vector<int> num_steps = {5, 10, 15, 20, 30, 50, 100};
    // unsigned int cap = 1000;
    // double precision = 1e-4;
    // double exact_value = 6.144;
    // char type = 'e';
    // plot_graph_performance(expiry, strike, cir_0, x_0, a, k , sigma, rho, r, type, num_steps, cap, precision, exact_value);

    return 0;
};
