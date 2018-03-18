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


int main(){
    // Initialize the random numbers Generator
    std::random_device rd;
    auto seed = rd();
    std::mt19937_64 generator(seed);

    unsigned int n = 10;

    // std::cout<<"Simulation of 10 normal distribution variables N(3, 1)"<<std::endl;
    // normal<std::mt19937_64>* variable = new normal<std::mt19937_64>(3, 1);
    // for(unsigned int i = 0; i<n; ++i){
    //     std::cout << (*variable)(generator) << std::endl;
    // };
    // std::cout << std::endl;
    //
    // std::cout << "Test of monte carlo on the gaussian N(3, 1), precision = 0.1" << std::endl;
    // monte_carlo<std::mt19937_64> mc(variable);
    // mc.set_precision(0.1);
    // mc.compute(generator);
    // mc.print();
    // std::cout << std::endl;
    //
    // std::cout << "Test of monte carlo and process on the integral of a brownian motion b/w 0 and 1."<< std::endl;
    // integral_brownian<std::mt19937_64>* ib = new integral_brownian<std::mt19937_64>(1000);
    // monte_carlo<std::mt19937_64> mc_2(ib);
    // mc_2.set_precision(0.1);
    // mc_2.compute(generator);
    // mc_2.print();
    //
    // std::cout << "Simulation of 10 variable normal_five_moments" << std::endl;
    // normal_seven_moments<std::mt19937_64> variable_5 = normal_seven_moments<std::mt19937_64>();
    // for(unsigned int i = 0; i<n; ++i){
    //     std::cout << variable_5(generator) << std::endl;
    // };
    // std::cout << std::endl;


    /// Plot CIR
    // double state_0 = 2;
    // double k = 1;
    // double a = 1;
    // double s = 1;
    // cir_o2<std::mt19937_64> c = cir_o2<std::mt19937_64>(state_0, k, a, s);
    // // cir_glasserman<std::mt19937_64> c = cir_glasserman<std::mt19937_64>(state_0, k, a, s);
    // c.set_time_end(1);
    // c.set_num_steps(10000);
    // c.write_to_plot(generator, 10, 0);
    // c.print_trajectory(0);

    // /// Plot brownian

    // brownian<std::mt19937_64> b = brownian<std::mt19937_64>();
    // b.set_num_steps(1000);
    // b.write_to_plot(generator, 20, 0);

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


    ///Plotting graphs to measure the performances of the different schemes
    double k = 0.5;
    double a = 0.02;
    double sigma = 0.4;
    double x_0 = 100.;
    double cir_0 = .04;
    double rho = -0.5;
    double r = 0.02;
    double strike = 100.;
    double maturity = 1.;
    char type = 'e';
    option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > * opt =
      new option<std::mt19937_64, heston<std::mt19937_64, cir_o2<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);

      option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > * opt_3 =
        new option<std::mt19937_64, heston<std::mt19937_64, cir_o3<std::mt19937_64> > > (maturity, strike, cir_0, x_0, a, k , sigma, rho, r, type);

    std::vector<int> num_steps = {5, 10, 15, 20, 30, 50, 100};
    int cap = 1e5;
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
      std::vector<double> temp = std::vector<double>();
      temp.push_back(mc.empirical_mean);
      temp.push_back(mc.empirical_mean - mc.ci_l_bound);

      //CIR_03
      opt_3->set_num_steps(num_steps.at(i));
      monte_carlo<std::mt19937_64> mc_3(opt_3);
      mc_3.set_precision(precision);
      mc_3.set_cap(cap);
      mc_3.compute(generator);
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

    return 0;
};
