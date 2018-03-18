#pragma once
#include "random_variable.h"

template<typename Generator> class normal : public random_variable<Generator>
{
public:
    normal();
    ~normal();
    normal(double m, double s);
    double operator()(Generator & gen);
private:
    double mean;
    double var;
    double std;
    std::normal_distribution<double> norm;
};

template <typename Generator> normal<Generator>::normal(){
    normal(0.0, 1.0);
};

template <typename Generator> normal<Generator>::normal(double m, double s){
    this->mean = m;
    this->var = std::pow(s, 2);
    this->std = s;
    this->norm= std::normal_distribution<double>(0.0,1.0);
};

template <typename Generator> double normal<Generator>::operator() (Generator & gen){
    this->realised = true;
    this->realisation = mean + std * norm(gen);
    return this->realisation;
};

template <typename Generator> normal<Generator>::~normal(){

};

const double x_1 = - std::sqrt(3);
const double x_2 = 0;
const double x_3 = std::sqrt(3);
const double p_1 = 1. / 6.;
const double p_2 = 5. / 6.;


template<typename Generator> class normal_five_moments : public random_variable<Generator>
{
public:
    normal_five_moments();
    ~normal_five_moments();

    double operator() (Generator & gen);
private:
    std::uniform_real_distribution<double> u;
};

template<typename Generator> normal_five_moments<Generator>::normal_five_moments(){
    this->u = std::uniform_real_distribution<double>(0, 1);
};

template<typename Generator> normal_five_moments<Generator>::~normal_five_moments(){
};

template<typename Generator> double normal_five_moments<Generator>::operator()(Generator & gen){
    double p = this->u(gen);
    if (p < p_1){
        return x_1;
    }
    else {
        if (p < p_2){
            return x_2;
        }
        else {
            return x_3;
        };
    };
};


const double y_1 = - std::sqrt(3. + std::sqrt(6.));
const double y_2 = - std::sqrt(3. - std::sqrt(6.));
const double y_3 = std::sqrt(3. - std::sqrt(6.));
const double y_4 = std::sqrt(3. + std::sqrt(6.));
const double q_1 = ( std::sqrt(6.) - 2. ) / ( 4. * std::sqrt(6.));
const double q_2 = 0.5;
const double q_3 = 1. - q_1;


template<typename Generator> class normal_seven_moments : public random_variable<Generator>
{
public:
    normal_seven_moments();
    ~normal_seven_moments();

    double operator() (Generator & gen);
private:
    std::uniform_real_distribution<double> u;
};

template<typename Generator> normal_seven_moments<Generator>::normal_seven_moments(){
    this->u = std::uniform_real_distribution<double>(0., 1.);
};

template<typename Generator> normal_seven_moments<Generator>::~normal_seven_moments(){
};

template<typename Generator> double normal_seven_moments<Generator>::operator()(Generator & gen){
    double q = (this->u)(gen);
    if (q < q_1){
        return y_1;
    }
    else {
        if (q < q_2){
            return y_2;
        }
        else {
            if (q < q_3){
                return y_3;
            }
            else {
                return y_4;
            };
        };
    };
};

const double r_1 = 1./3.;
const double r_2 = 2./3.;

template<typename Generator> class zeta : public random_variable<Generator>
{
public:
  zeta();
  ~zeta();

  double operator() (Generator & gen);
private:
  std::uniform_real_distribution<double> u;
};

template<typename Generator> zeta<Generator>::zeta(){
  this->u = std::uniform_real_distribution<double>(0., 1.);
};

template<typename Generator> zeta<Generator>::~zeta(){
};

template<typename Generator> double zeta<Generator>::operator() (Generator & gen){
  double q = (this->u)(gen);
  if (q < r_1){
    return 1.;
  }
  else{
    if (q < r_2){
      return 2.;
    }
    else{
      return 3.;
    }
  }
};


template<typename Generator> class rademacher : public random_variable<Generator>
{
public:
  rademacher();
  ~rademacher();

  double operator() (Generator & gen);
private:
  std::uniform_real_distribution<double> u;
};

template<typename Generator> rademacher<Generator>::rademacher(){
  this->u = std::uniform_real_distribution<double>(0., 1.);
};

template<typename Generator> rademacher<Generator>::~rademacher(){
};

template<typename Generator> double rademacher<Generator>::operator() (Generator & gen){
  double q = (this->u)(gen);
  if (q < 0.5){
    return 1.;
  }
  else{
    return -1.;
  }
};
