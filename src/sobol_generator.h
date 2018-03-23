#pragma once
#include <random>
#include <iostream>
#include <math.h>
#include <utility>

#include <gsl/gsl_qrng.h>

class sobol_generator
{
public:
    sobol_generator();
    ~sobol_generator();

    double operator()();
    double min();
    double max();
private:
    gsl_qrng* q;
    double* v;
};

sobol_generator::sobol_generator(){
    this->q = gsl_qrng_alloc (gsl_qrng_sobol, 2);
    this->v = new double[1];
};

sobol_generator::~sobol_generator(){
    gsl_qrng_free(this->q);
    delete this->v;
};

double sobol_generator::operator()(){
    gsl_qrng_get (q, v);
    return v[0];
};

double sobol_generator::min(){
    return 0.;
};

double sobol_generator::max(){
    return 1.;
};
