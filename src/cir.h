#pragma once
#include "random_variable.h"
#include "normal.h"
#include "process.h"

template <typename Generator> class cir : public process<Generator>
{
    /// Generation of a class with a CIR process
    /// dX = (a - kX) dt + sigma sqrt(X) dW
    /// X(0) = x_0
public:
    /// Constructors & destructors
    cir();
    cir(double state_0, double k, double a, double s);
    ~cir();
    // cir(cir<Generator> & c
    virtual std::vector<double> next_step(Generator & gen, std::vector<double> state, double t) = 0;
    /// Getters and Setters
    double get_a();
    double get_sigma();
    void set_a(double a);
    void set_sigma(double s);
protected:
    double a;
    double sigma;
    double k;
};

template <typename Generator> cir<Generator>::cir(): cir(1., 0., 1., 1.) {
};

template <typename Generator> cir<Generator>::~cir(){
};

template <typename Generator> cir<Generator>::cir(double state_0, double k, double a, double s): process<Generator>()
{
    this->dimension = 1.;
    this->state_0 = std::vector<double> (this->dimension, state_0);
    this->a = a;
    this->k = k;
    this->sigma = s;
    this->set_state_0(this->state_0);
    this->state = this->state_0;
};

template <typename Generator> double cir<Generator>::get_a(){
    return this->a;
};

template <typename Generator> double cir<Generator>::get_sigma(){
    return this->sigma;
};

template <typename Generator> void cir<Generator>::set_a(double a){
    this->a = a;
};

template <typename Generator> void cir<Generator>::set_sigma(double s){
    this->sigma = s;
};

//////////
/// CIR_O2
//////////
template <typename Generator> class cir_o2 : public cir<Generator>
{
    /// Generation of a class with a CIR process
    /// dX = (a - kX) dt + sigma sqrt(X) dW
    /// X(0) = x_0
public:
    /// Constructors & destructors
    cir_o2();
    cir_o2(double state_0, double k, double a, double s);
    ~cir_o2();
    // cir(cir<Generator> & c);

    std::vector<double> next_step(Generator & gen, std::vector<double> state, double t);
    /// Getters and Setters
protected:
    std::uniform_real_distribution<double> u;
    normal_five_moments<Generator> y;
};

double psi_f(double k, double t){
    if (k != 0){
        return ((1 - std::exp(-k * t)) / k);
    }
    else {
        return t;
    };
};

double phi_f(double x, double t, double w, double sigma, double a, double k){
    double expo = std::exp(- k * t / 2);
    double psi = psi_f(k, t/2);
    double aux = a - std::pow(sigma, 2)/4;
    double result = expo;
    result *= std::pow( std::sqrt(aux * psi + expo * x) + sigma * w / 2 , 2);
    result += aux * psi;
    return result;
};


double k_2_f(double t, double sigma, double a, double k){
    double expo = std::exp(+ k * t / 2);
    double psi = psi_f(k, t/2);
    double aux = std::pow(sigma, 2)/4 - a;
    if (aux > 0){
        double result = expo;
        result *= aux * psi + std::pow( std::sqrt(expo * aux * psi) + sigma / 2 * std::sqrt(3 * t), 2);
        return result;
    }
    else {
        return 0;
    };
};

double u_tilde_1_f(double t, double x, double a, double k){
    double psi = psi_f(k, t);
    return x * std::exp(-k * t) + a * psi;
};

double u_tilde_2_f(double t, double x, double sigma, double a, double k){
    double psi = psi_f(k, t);
    double expo = std::exp(-k*t);
    double result = std::pow( u_tilde_1_f(t, x, a, k), 2) + std::pow(sigma, 2) * psi * (a * psi / 2 + x * expo);
    return result;
};

double pi_f(double u_1, double u_2){
    double delta = 1 - std::pow(u_1, 2) / u_2;
    double result = ( 1. - std::sqrt(delta) ) / 2.;
    return result;
};


template<typename Generator> cir_o2<Generator>::cir_o2(): cir<Generator> ()
{
    this->u = std::uniform_real_distribution<double>(0., 1.);
};
template<typename Generator> cir_o2<Generator>::cir_o2(double state_0, double k, double a, double s): cir<Generator> (state_0, k, a, s)
{
    this->u = std::uniform_real_distribution<double>(0., 1.);
    this->y = normal_five_moments<Generator>();
};


template<typename Generator> cir_o2<Generator>::~cir_o2(){
};

template<typename Generator> std::vector<double> cir_o2<Generator>::next_step(Generator & gen, std::vector<double> state, double t) {
    double k_2 = k_2_f(t, this->sigma, this->a, this->k);
    double x = state.at(0);
    double tt = this->time_step;
    if( x >= k_2){
        double normal_test = (this->y)(gen);
        double w = std::sqrt(tt) * normal_test;
        return std::vector<double> (this->dimension, phi_f(x, tt, w, this->sigma, this->a, this->k));
    }
    else {
        double u_1 = u_tilde_1_f(tt, x, this->a, this->k);
        double u_2 = u_tilde_2_f(tt, x, this->sigma, this->a, this->k);
        double pi = pi_f(u_1, u_2);
        double uniform = (this->u)(gen);
        if(uniform < pi){
            return std::vector<double> (this->dimension, u_1 / (2. * pi) );
        }
        else {
            return std::vector<double> (this->dimension , u_1 / (2. * (1. - pi) ));
        };
    };
};

//////////////////
/// CIR_O3
//////////////////
template <typename Generator> class cir_o3 : public cir<Generator>
{
    /// Generation of a class with a CIR process
    /// dX = (a - kX) dt + sigma sqrt(X) dW
    /// X(0) = x_0
public:
    /// Constructors & destructors
    cir_o3();
    cir_o3(double state_0, double k, double a, double s);
    ~cir_o3();
    // cir(cir<Generator> & c);

    std::vector<double> next_step(Generator & gen, std::vector<double> state, double t);
    /// Getters and Setters
protected:
    std::uniform_real_distribution<double> u;
    normal_seven_moments<Generator> seven;
    zeta<Generator> zet;
    rademacher<Generator> rad;
};

double X_0_f(double x, double t, double sigma, double a, double k){
  double aux = a - std::pow(sigma, 2) / 4.;
  double psi = psi_f(-k, t);
  return x + aux * psi;
};

double X_1_f(double x, double t, double sigma, double k, double y){
  double psi = psi_f(-k, t);
  double temp = std::sqrt(x) + sigma * std::sqrt(psi) * y / 2.;
  double result = 0.;
  if (temp > 0){
    result = std::pow(temp, 2);
  }
  return result;
};

double X_t_f(double x, double t, double sigma, double a, double k, double eps){
  double aux = a - std::pow(sigma, 2) / 4.;
  double psi = psi_f(-k, t);
  if (aux < 0){
    aux = -aux;
  }
  double result = x;
  result += sigma / std::sqrt(2.) * std::sqrt(aux) * eps * psi;
  return result;
};

double k_3_f(double t, double sigma, double a, double k){
  double aux = a - std::pow(sigma, 2) / 4.;
  double psi = psi_f(-k, t);
  double num = std::sqrt(3. + std::sqrt(6.));
  double temp;
  if (std::pow(sigma, 2) <= 4 * a / 3.){
    temp = sigma / std::sqrt(2.) * std::sqrt(aux);
  }
  else{
    if (aux > 0){
      temp = sigma / std::sqrt(2.) * num;
      temp += std::sqrt(- aux + sigma / std::sqrt(2.) * std::sqrt(aux));
      temp = std::pow(temp, 2);
    }
    else{
      temp = sigma / std::sqrt(2.) * num;
      temp += std::sqrt(sigma / std::sqrt(2.) * std::sqrt(- aux));
      temp = std::pow(temp, 2);
      temp += (- aux);
    }
  }
  return psi * temp;
};

double u_tilde_3_f(double t, double x, double sigma, double a, double k){
  double u_1 = u_tilde_1_f(t, x, a, k);
  double u_2 = u_tilde_2_f(t, x, sigma, a, k);
  double psi = psi_f(k, t);
  double expo = std::exp(- k * t);

  double factor = 3 * x * expo + a * psi;
  factor *= psi * (a + std::pow(sigma, 2) / 2.);
  factor += 2 * std::pow(x, 2) * std::pow(expo, 2);
  factor *= std::pow(sigma, 2) * psi;

  double result = u_1 * u_2;
  result += factor;
  return result;
};

template<typename Generator> cir_o3<Generator>::cir_o3(): cir<Generator> ()
{
    this->u = std::uniform_real_distribution<double>(0., 1.);
};
template<typename Generator> cir_o3<Generator>::cir_o3(double state_0, double k, double a, double s): cir<Generator> (state_0, k, a, s)
{
    this->u = std::uniform_real_distribution<double>(0., 1.);
    this->seven = normal_seven_moments<Generator>();
    this->zet = zeta<Generator>();
    this->rad = rademacher<Generator>();
};

template<typename Generator> cir_o3<Generator>::~cir_o3(){
};

template<typename Generator> std::vector<double> cir_o3<Generator>::next_step(Generator & gen, std::vector<double> state, double t){
  double result = 0.;
  double zeta = (this->zet)(gen);
  double y = (this->seven)(gen);
  double eps = (this->rad)(gen);
  double k_3 = k_3_f(t, this->sigma, this->a, this->k);
  double x = state.at(0);
  double tt = this->time_step;
  double aux = 4 * this->a - std::pow(this->sigma, 2);
  if (x >= k_3){
    result = x;
    if (zeta == 1.){
      if (aux >= 0){
        result = X_1_f(result, tt, this->sigma, this->k, y);
        result = X_0_f(result, tt, this->sigma, this->a, this->k);
        result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
      }
      else{
        result = X_0_f(result, tt, this->sigma, this->a, this->k);
        result = X_1_f(result, tt, this->sigma, this->k, y);
        result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
      }
    }
    else{
      if (zeta == 2.){
        if (aux >= 0){
          result = X_1_f(result, tt, this->sigma, this->k, y);
          result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
          result = X_0_f(result, tt, this->sigma, this->a, this->k);
        }
        else{
          result = X_0_f(result, tt, this->sigma, this->a, this->k);
          result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
          result = X_1_f(result, tt, this->sigma, this->k, y);
        }
      }
      else{
        if (aux >= 0){
          result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
          result = X_1_f(result, tt, this->sigma, this->k, y);
          result = X_0_f(result, tt, this->sigma, this->a, this->k);
        }
        else{
          result = X_t_f(result, tt, this->sigma, this->a, this->k, eps);
          result = X_0_f(result, tt, this->sigma, this->a, this->k);
          result = X_1_f(result, tt, this->sigma, this->k, y);
        }
      }
    }
    result *= std::exp(- this->k * tt);
  }
  else{
    double u_1 = u_tilde_1_f(tt, x, this->a, this->k);
    double u_2 = u_tilde_2_f(tt, x, this->sigma, this->a, this->k);
    double u_3 = u_tilde_3_f(tt, x, this->sigma, this->a, this->k);
    double s = (u_3 - u_1 * u_2) / (u_2 - std::pow(u_1, 2));
    double p = (u_1 * u_3 - std::pow(u_2, 2)) / (u_2 - std::pow(u_1, 2));
    double delta = std::sqrt(std::pow(s, 2) - 4 * p);
    double pi = (u_1 - (s - delta) / 2.) / delta;
    double uniform = (this->u)(gen);
    if (uniform < pi){
      result = (s + delta) / 2.;
    }
    else{
      result = (s - delta) / 2.;
    }
  }

  return std::vector<double> (this->dimension, result);
};


//////////////////
/// CIR_glasserman
//////////////////
template <typename Generator> class cir_glasserman : public cir<Generator>
{
    /// Generation of a class with a CIR process
    /// dX = (a - kX) dt + sigma sqrt(X) dW
    /// X(0) = x_0
public:
    /// Constructors & destructors
    cir_glasserman();
    cir_glasserman(double state_0, double k, double a, double s);
    ~cir_glasserman();
    // cir(cir<Generator> & c);

    std::vector<double> next_step(Generator & gen, std::vector<double> state, double t);
    /// Getters and Setters
protected:
    std::poisson_distribution<int> n_poisson;
    std::chi_squared_distribution<double> x_chi;
    std::normal_distribution<double> z_norm;
    double d;
};

template<typename Generator> cir_glasserman<Generator>::cir_glasserman(): cir_glasserman<Generator> (1, 1, 1, 1)
{

};

template<typename Generator> cir_glasserman<Generator>::cir_glasserman(double state_0, double k, double a, double s): cir<Generator> (state_0, k, a, s)
{
    if(this->k == 0){
        std::cout << "We do not accept k = 0 for Glasserman's simulations" << std::endl;
        return ;
    };
    this->d = 4. * this->a / std::pow(this->sigma, 2);
    this->n_poisson = std::poisson_distribution<int>(1.);
    this->x_chi = std::chi_squared_distribution<double> (d-1);
    this->z_norm = std::normal_distribution<double> (0., 1.);
};

template<typename Generator> cir_glasserman<Generator>::~cir_glasserman(){

};

template<typename Generator> std::vector<double> cir_glasserman<Generator>::next_step(Generator & gen, std::vector<double> state, double t) {
    double x = state.at(0);
    double c  = std::pow(this->sigma, 2) * (1. - std::exp(-this->k * this->time_step) ) / 4. / this->k;
    double lam = x * std::exp(- this->k * this->time_step) / c;
    if(this->d <= 1){
        n_poisson = std::poisson_distribution<int>(lam/2.);
        int poisson = n_poisson(gen);
        x_chi = std::chi_squared_distribution<double>(d + 2. * poisson);
        double chi = x_chi(gen);
        return std::vector<double> (this->dimension, c * chi);
    }
    else {
        double chi = x_chi(gen);
        double norm = z_norm(gen);
        return std::vector<double> (this->dimension, c * ( std::pow( (norm + std::sqrt(lam)), 2 ) + chi ) );
    };
};
