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
    unsigned int num_steps = 100;
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


int main(){
    /// Try heston with Glasserman
    // double k = 0.5;
    // double a = 0.02;
    // double sigma = 0.2;
    // double x_0 = 100;
    // double cir_0 = .04;
    // double rho = -0.3;
    // double r = 0.02;
    // heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > h = heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    // h.set_num_steps(100);
    // h(generator);
    // h.write_to_plot(generator, 10, 2);
    //h.print_trajectory(2);

    //Try Option with Glasserman
    // double strike = 100.;
    // double maturity = 1.;
    // char type = 'a';
    // option<std::mt19937_64, heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > > * opt =
    //   new option<std::mt19937_64, heston<std::mt19937_64, cir_glasserman<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // for(unsigned int i = 0; i<n; ++i){
    //     std::cout << (*opt)(generator) << std::endl;
    // };
    // std::cout << std::endl;
    //
    // std::cout << "Test of Monte-Carlo on an European option on price" << std::endl;
    // monte_carlo<std::mt19937_64> mc(opt);
    // mc.set_precision(0.1);
    // mc.set_cap(10000);
    // mc.compute(generator);
    // mc.print();
    // std::cout << std::endl;

    /// Try heston with CIR_O2
    // The parameters are chosen to reproduce the graphs of Alfonsi p.26
    // double k = 0.5;
    // double a = 0.02;
    // double sigma = 0.4;
    // double x_0 = 100;
    // double cir_0 = .04;
    // double rho = -0.5;
    // double r = 0.02;
    //heston<std::mt19937_64, cir_o2<std::mt19937_64> > h = heston<std::mt19937_64, cir_o2<std::mt19937_64> > (cir_0, x_0, a, k , sigma, rho, r);
    //h.set_num_steps(10);
    //h(generator);
    //h.write_to_plot(generator, 10, 2);
    //h.print_trajectory(2);

    //Try Option with CIR_O2
    // double k = 0.5;
    // double a = 0.02;
    // double sigma = 0.4;
    // double x_0 = 100;
    // double cir_0 = .04;
    // double rho = -0.5;
    // double r = 0.02;
    // double strike = 100.;
    // double maturity = 1.;
    // char type = 'e';
    // option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt =
    //   new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // opt->set_num_steps(10);
    //
    // std::cout << "Test of Monte-Carlo on an European option (put) on price" << std::endl;
    // monte_carlo<std::mt19937_64> mc(opt);
    // mc.set_precision(0.01);
    // mc.set_cap(1000000);
    // mc.compute(generator);
    // mc.print();
    // std::cout << std::endl;

    //Try Option with CIR_O3
    // double k = 0.5;
    // double a = 0.02;
    // double sigma = 0.4;
    // double x_0 = 100;
    // double cir_0 = .04;
    // double rho = -0.5;
    // double r = 0.02;
    // double strike = 100.;
    // double maturity = 1.;
    // char type = 'e';
    // option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > * opt =
    //   new option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // opt->set_num_steps(20);
    //
    // std::cout << "Test of Monte-Carlo on an European option (put) on price" << std::endl;
    // monte_carlo<std::mt19937_64> mc(opt);
    // mc.set_precision(0.1);
    // mc.set_cap(100000);
    // mc.compute(generator);
    // mc.print();
    // std::cout << std::endl;

    //Try Closed formula for first order moment of log spot (log_spot_one)
    // double k = 0.5;
    // double a = 0.02;
    // double sigma = 2.0;
    // double x_0 = 100.;
    // double cir_0 = .04;
    // double rho = -1.;
    // double r = 0.02;
    // double strike = 100.;
    // double maturity = 1.;
    // char type = 'e';
    // option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt =
    //   new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // opt->set_num_steps(20);
    //
    // double theo = opt->theoretical_log_spot_one();
    //
    // std::cout<<"Test for the closed formula of first order moment of log spot" <<std::endl;
    // std::cout << "Theoretical Value" << "\t" << theo << std::endl;
    // monte_carlo<std::mt19937_64> mc(opt);
    // mc.set_precision(0.01);
    // mc.set_cap(10000);
    // mc.compute(generator);
    // mc.print();
    // std::cout << std::endl;
    //
    //



    // option<sobol_generator, heston<sobol_generator, cir_o2<sobol_generator> > > * opt =
    //   new option<sobol_generator, heston<sobol_generator, cir_o2<sobol_generator> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // // // option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > * opt_3 =
    // // //     new option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);
    //
    // unsigned int cap = 3*1e5;
    // unsigned num_steps = 100;
    // double precision = 1e-5;

    // opt->set_num_steps(num_steps);
    // monte_carlo<sobol_generator> mc(opt);
    // mc.set_precision(precision);
    // mc.set_cap(cap);
    // std::cout << "Monte Carlo without control variance on CIR_O2 Heston Asian" << std::endl;
    // mc.compute(sobol);
    // mc.print();

    //
    // normal<sobol_generator>* norm_sobol = new normal<sobol_generator>(3, 5);
    // monte_carlo<sobol_generator> mc_norm(norm_sobol);
    // mc_norm.set_precision(precision);
    // mc_norm.set_cap(cap);
    // std::cout << "Monte Carlo without control variance on N(3, 5)" << std::endl;
    // mc_norm.compute(sobol);
    // mc_norm.print();
    //
    //
    // normal_five_moments<sobol_generator>* norm_5 = new normal_five_moments<sobol_generator>();
    // monte_carlo<sobol_generator> mc_5(norm_5);
    // mc_5.set_precision(precision);
    // mc_5.set_cap(cap);
    // std::cout << "Monte Carlo without control variance on Normal five moments / Should work" << std::endl;
    // mc_5.compute(sobol);
    // mc_5.print();
    //

    // std::cout << std::endl;
    // std::cout << "Monte Carlo with control variance" <<std::endl;
    // mc.compute_control_variate(generator);
    // mc.print();

/*
    std::vector<int> num_steps = {5, 10, 15, 20, 30, 50, 100};
    int cap = 1*1e3;
    double precision = 1e-4;
    double exact_value = 6.144;
    std::vector<std::vector<double> > data = std::vector<std::vector<double> > ();

    for (unsigned int i = 0; i < num_steps.size(); ++i){
      //CIR_O2
      opt->set_num_steps(num_steps.at(i));
      monte_carlo<std::mt19937_64> mc(opt);
      mc.set_precision(precision);
      mc.set_cap(cap);
      mc.compute(generator);
    //   mc.compute_control_variate(generator);
      std::vector<double> temp = std::vector<double>();
      temp.push_back(mc.empirical_mean);
      temp.push_back(mc.empirical_mean - mc.ci_l_bound);

      //CIR_03
      opt_3->set_num_steps(num_steps.at(i));
      monte_carlo<std::mt19937_64> mc_3(opt_3);
      mc_3.set_precision(precision);
      mc_3.set_cap(cap);

      mc_3.compute(generator);
    //   mc_3.compute_control_variate(generator);

      temp.push_back(mc_3.empirical_mean);
      temp.push_back(mc_3.empirical_mean - mc_3.ci_l_bound);

      data.push_back(temp);
    }
    std::cout << "We have computed the trajectories. " << std::endl;

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
*/


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

    /// Plotting som examples of trajectories.
    /// The last parameter is the coordinate : 0 = CIR, 2 = Heston
    // plot_heston_glasserman(cir_0, x_0, a, k, sigma, rho, r, num_steps, 2);
    // plot_heston_o2(cir_0, x_0, a, k, sigma, rho, r, num_steps, 0);
    // plot_heston_o3(cir_0, x_0, a, k, sigma, rho, r, num_steps, 2);


    /// Testing the control variates
    // compare_reduction_variance_glasserman(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'e', cap);
    // compare_reduction_variance_glasserman(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'a', cap);

    // compare_reduction_variance_o2(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'e', cap);
    // compare_reduction_variance_o3(cir_0, x_0, a, k, sigma, rho, r, strike, expiry, num_steps, 'e', cap);

    /// One shot to compare with cuda code
    // test_for_cuda();

    /// Test for Sobol
    // test_mc_sobol(cir_0, x_0, a, k, sigma, rho, r, strike,  expiry, num_steps, cap);

    return 0;
};
