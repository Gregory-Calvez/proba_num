#pragma once
#include "random_variable.h"
#include "process.h"
#include "cir.h"

template<typename Generator, typename Heston> class option : public random_variable<Generator>
{
public:
  //Constructors and destructors
  option();
  ~option();
  option(double maturity, double strike, double cir_0, double x_0,
    double a, double k, double sigma, double rho, double r, char type);

  double payoff(std::vector<std::vector<double>> trajectory);
  double test_log_spot_one(std::vector<std::vector<double>> trajectory);
  double operator()(Generator & gen);
  void set_num_steps(unsigned int n);
  double theoretical_log_spot_one();


private:
  //Option parameters
  double maturity;
  double strike;
  char type;

  //Heston parameters
  double r;
  double a;
  double k;
  double sigma;
  double rho;
  Heston heston;
};

template<typename Generator, typename Heston> option<Generator, Heston>::~option()
{
};

template<typename Generator, typename Heston> option<Generator, Heston>::option() : option<Generator, Heston>(1.0, 100.,1, 100, 1, 1, 1, 0, 0 )
{
};

template<typename Generator, typename Heston> option<Generator, Heston>::option(double maturity, double strike, double cir_0, double x_0,
  double a, double k, double sigma, double rho, double r, char type) : random_variable<Generator>(), heston(cir_0, x_0, a, k, sigma, rho, r), maturity(maturity),
  strike(strike), type(type), r(r), a(a), k(k), sigma(sigma), rho(rho)
{
};

template<typename Generator, typename Heston> void option<Generator, Heston>::set_num_steps(unsigned int n)
{
  (this->heston).set_num_steps(n);
};

template<typename Generator, typename Heston> double option<Generator, Heston>::payoff(std::vector<std::vector<double>> trajectory){
  int indice = 0;
  if (type == 'e'){
    indice = 2;
  }
  else {
    if (type == 'a'){
      indice = 3;
    }
    else {
      std::cout << "Option type not recognised" << std::endl;
    }
  }
  double r = this->r;
  double T = this->maturity;
  double K = this->strike;
  int num_steps = (this->heston).get_num_steps();
  std::vector<double> final_state = trajectory.at(num_steps);
  double discount = std::exp(- r * T);
  double X = final_state.at(indice);
  double result = 0.0;
  if (X < K){
    result += K - X;
  }
  result *= discount;
  return result;
};

//Real version
template<typename Generator, typename Heston> double option<Generator, Heston>::operator()(Generator & gen){
  std::vector<std::vector<double>> trajectory = (this->heston)(gen);
  this->realised = true;
  this->realisation = this->payoff(trajectory);
  return this->realisation;
};

//Version to test log_spot_one
// template<typename Generator, typename Heston> double option<Generator, Heston>::operator()(Generator & gen){
//   std::vector<std::vector<double>> trajectory = (this->heston)(gen);
//   this->realised = true;
//   this->realisation = this->test_log_spot_one(trajectory);
//   return this->realisation;
// };

template<typename Generator, typename Heston> double option<Generator, Heston>::test_log_spot_one(std::vector<std::vector<double>> trajectory){
  double result;
  double r = this->r;
  double T = this->maturity;
  int num_steps = (this->heston).get_num_steps();
  std::vector<double> final_state = trajectory.at(num_steps);
  double S = final_state.at(2);
  result += std::log(S) - r * T;
  return result;
};

template<typename Generator, typename Heston> double option<Generator, Heston>::theoretical_log_spot_one(){
  double T = this->maturity;
  double theo = (this->heston).log_spot_one(T);
  return theo;
};
