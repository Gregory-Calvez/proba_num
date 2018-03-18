#pragma once
#include "process.h"
#include "cir.h"

template<typename Generator, typename Cir> class heston : public process<Generator>
{
public:
    /// Constructors & destructors
    heston();
    ~heston();
    heston(double cir_0, double x_0, double a, double k, double sigma, double rho, double r);
    // heston(heston<Generator & h);

    /// Getters and setters

    /// Functions
    std::vector<double> hw(Generator & gen, std::vector<double> state, double t);
    std::vector<double> hz(Generator & gen, std::vector<double> state, double t);
    std::vector<double> next_step(Generator & gen, std::vector<double> state, double t);
    double log_spot_one(double t);

private:
    // X_1, .., X_4 pour tout temps.
    // We still keep only the stock price, which is of interest
    // in the process member.
    double x_0;
    double cir_0;
    double r;
    double a;
    double k;
    double sigma;
    double rho;
    std::bernoulli_distribution bernoulli;
    std::normal_distribution<double> norm;
    Cir cir;
};


template<typename Generator, typename Cir> heston<Generator, Cir>::heston() : heston<Generator, Cir>(1, 100, 1, 1, 1, 0, 0)
{
};

template<typename Generator, typename Cir> heston<Generator, Cir>::heston(double cir_0, double x_0, double a, double k, double sigma, double rho, double r)
    : process<Generator>(), cir(cir_0, k, a, sigma), x_0(x_0), cir_0(cir_0), r(r), a(a), k(k), sigma(sigma), rho(rho)
{
    // TODO constructor w/ dimension for process
    this->dimension = 4;
    this->state_0 = std::vector<double> ({cir_0, 0, x_0, 0});
    this->state = std::vector<double> ({cir_0, 0, x_0, 0});
    norm = std::normal_distribution<double>(0, 1);
    bernoulli = std::bernoulli_distribution(0.5);
};

template<typename Generator, typename Cir> heston<Generator, Cir>::~heston(){
};

template<typename Generator, typename Cir> std::vector<double> heston<Generator, Cir>::hw(Generator & gen, std::vector<double> state, double t){
    double tt = this->time_step;
    double x_1 = state.at(0);
    double x_2 = state.at(1);
    double x_3 = state.at(2);
    double x_4 = state.at(3);
    double d_x_1 = - x_1;
    std::vector<double> cir_update = cir.next_step(gen, std::vector<double> (1, state.at(0)), tt);
    cir.set_state(cir_update);
    d_x_1 = d_x_1 + cir_update.at(0);
    x_2 = x_2 + (x_1 + 0.5 * d_x_1) * tt;
    x_4 = x_4 + 0.5 * x_3 * tt;
    x_3 = x_3 * std::exp( (this->r - this->rho * this->a / this->sigma) * tt + this->rho * d_x_1 / this->sigma + (this->rho * this-> k / this->sigma - 0.5) * (x_1 + 0.5 * d_x_1) * tt );
    x_4 = x_4 + 0.5 * x_3 * tt;
    x_1 = x_1 + d_x_1;
    this->state = std::vector<double>({x_1, x_2, x_3, x_4});
    return this->state;
};
template<typename Generator, typename Cir> std::vector<double> heston<Generator, Cir>::hz(Generator & gen, std::vector<double> state, double t){
    double tt = this->time_step;
    double x_1 = state.at(0);
    double x_2 = state.at(1);
    double x_3 = state.at(2);
    double x_4 = state.at(3);
    double normal_realisation = norm(gen);
    x_3 = x_3 * std::exp(normal_realisation * std::sqrt( (1-std::pow(rho, 2)) * x_1 * tt ));
    this->state = std::vector<double>({x_1, x_2, x_3, x_4});
    return this->state;
};
template<typename Generator, typename Cir> std::vector<double> heston<Generator, Cir>::next_step(Generator & gen, std::vector<double> state, double t)
{
    //std::cout << this->state.at(0) << " " << this->state.at(1) << " " << this->state.at(2) << " " << this->state.at(3) << std::endl;
    bool bernoulli_realisation = bernoulli(gen);
    if(bernoulli_realisation){
        this->hz(gen, this->state, t);
        this->hw(gen, this->state, t);
    }
    else {
        this->hw(gen, this->state, t);
        this->hz(gen, this->state, t);
    }
    return this->state;
};

template<typename Generator, typename Cir> double heston<Generator, Cir>::log_spot_one(double t){
  double x_0 = this->x_0;
  double cir_0 = this->cir_0;
  double a = this -> a;
  double k = this->k;
  double sigma = this->sigma;
  double rho = this->rho;

  double sinus = std::sinh(k * t / 2.);
  double cosinus = std::cosh(k * t / 2.);
  double expo = std::exp(- k * t / 2.);
  double p_prime_zero = std::pow(sigma, 2) / (2. * k) - sigma * rho;

  double result = std::log(x_0);
  double temp = t / 2. * p_prime_zero - sigma * rho / k - p_prime_zero / k;

  temp *= sinus;
  temp += t / 2. * p_prime_zero * cosinus;
  temp *= - expo;
  temp += - sigma * rho * t / 2.;
  temp *= 2 * a / std::pow(sigma, 2);

  result += temp;
  result += -cir_0 * expo * sinus / k;

  return result;
};
