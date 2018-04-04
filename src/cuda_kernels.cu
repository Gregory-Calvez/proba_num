#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>

#include <time.h>

const int num_blocks = 1024;
const int num_threads = 256;
const int num_iterations = 10;


__global__ void setup_states(curandState* states){
    int id = threadIdx.x + num_threads * blockIdx.x;
    // Initialisation states
    curand_init(0, id, 0, &states[id]);
}

__global__ void put_euro_heston_cir_o2(curandState* states, float* results_sum, float* results_sum_squared, float x_0, float cir_0, float r, float a, float k, float sigma, float rho, float t, float strike, unsigned int num_steps, char type)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    // Saving the state in the GPU memory to be more efficient
    curandState localState = states[id];
    // Shared memory for the Monte Carlo
    __shared__ float partial_sums[num_threads];
    __shared__ float partial_sums_squared[num_threads];

    // Some values that we will not have to compute twice / maybe we should give this in parameters.
    float x_1;
    float x_2;
    float x_3;
    float x_4;
    float dx_1;
    float dt = t / num_steps;
    float u_1; // u will be a uniform variable
    float u_2; //
    float y;
    float2 n; // n will be a pair of normal variable
    float u_tilde_1;
    float u_tilde_2;
    float pi;
    float value_option;

    float aux_hz_1 = (1 - rho*rho) * dt;
    float aux_hw_1 = (r - rho * a / sigma) * dt;
    float aux_hw_2 = rho / sigma;
    float aux_hw_3 = (rho*k/sigma - 0.5) * dt;
    float srqt3 = sqrtf(3);
    float aux_phi_1;

    float aux_sigma2sur4moinsa = sigma*sigma / 4 - a;
    float psi_k;
    if (k == 0){
        psi_k = dt / 2;
    }
    else {
        psi_k = ( 1 - expf(-k * dt / 2) ) / k;
    };
    float k_2;
    float expo = expf(k * dt / 2);
    float expo_2 = expf(-k*dt);
    float aux_k_2_1 = sqrtf(expo * aux_sigma2sur4moinsa * psi_k) + sigma / 2 * sqrt(3*t);
    if (aux_sigma2sur4moinsa > 0){
        k_2 = expo * (aux_sigma2sur4moinsa * psi_k + aux_k_2_1 * aux_k_2_1);
    }
    else {
        k_2 = 0;
    };

    float psi_k_2;
    if (k == 0){
        psi_k_2 = dt;
    }
    else {
        psi_k_2 = ( 1 - expf(-k * dt) ) / k;
    };
    float aux_u_tilde_1 = a*psi_k_2;

    // Main loop for the Monte Carlo
    for(unsigned int i = 0; i < num_iterations; ++i){
        // Initialization of the shared memory at the begining of the MC
        if(i == 0){
            partial_sums[threadIdx.x] = 0;
            partial_sums_squared[threadIdx.x] = 0;
        };
        x_1 = cir_0;  // Vol process
        x_2 = 0;  // Integration of the vol process
        x_3 = x_0;  // Stock process
        x_4 = 0;  // Integration of the stock process

        // CIR_O2 and Heston
        for (unsigned int k = 0; k < num_steps; ++k){
            // printf("Value of the Heston %f  and the CIR %f at the step %d. \n", x_3, x_1, k);
            u_1 = curand_uniform(&localState);
            n = curand_normal2(&localState); // It is not optimal, we simulate two uniform for 1 normal
            if (u_1 < 0.5) {
                // HZ
                x_3 = x_3 * expf(sqrtf(x_1 * aux_hz_1)*n.x);
                // HW
                dx_1 = - x_1;
/////////////// CIR_O2
                if(x_1 >= k_2){
                    u_2 = curand_uniform(&localState);
                    if(u_2 < 1./6.){
                        y = - srqt3;
                    }
                    else {
                        if (u_2 < 5./6.){
                            y = 0;
                        }
                        else {
                            y = srqt3;
                        };
                    };
                    aux_phi_1 = sqrtf( - aux_sigma2sur4moinsa * psi_k + x_1 /expo ) + sigma / 2 * sqrtf(dt) * y;
                    x_1 = 1 / expo * aux_phi_1 * aux_phi_1 - aux_sigma2sur4moinsa * psi_k;
                }
                else {
                    u_tilde_1 = x_1 * expo_2 + aux_u_tilde_1;
                    u_tilde_2 = u_tilde_1 * u_tilde_1 + sigma * sigma * psi_k_2 * (a * psi_k_2 / 2.0f + x_1 * expo_2);
                    pi = 0.5f * (1 - sqrtf(1 - u_tilde_1 * u_tilde_1 / u_tilde_2) );
                    u_2 = curand_uniform(&localState);
                    if (u_2 < pi){
                        x_1 = u_tilde_1 / 2.0f / pi;
                    }
                    else {
                        x_1 = u_tilde_1 / 2.0f / (1.0f - pi);
                    };
                };
/////////////// Fin CIR_O2

                dx_1 += x_1;
                x_2 += (x_1 - 0.5*dx_1) * dt;
                x_4 += 0.5*x_3*dt;
                x_3 = x_3 * expf( aux_hw_1 + aux_hw_2 * dx_1 + aux_hw_3 * (x_1 - 0.5 * dx_1) );
                x_4 += 0.5*x_3*dt;
            }
            else {
/////////////// HW
                dx_1 = - x_1;
/////////////// CIR_O2
                if(x_1 >= k_2){
                    u_2 = curand_uniform(&localState);
                    if(u_2 < 1./6.){
                        y = - srqt3;
                    }
                    else {
                        if (u_2 < 5./6.){
                            y = 0;
                        }
                        else {
                            y = srqt3;
                        };
                    };
                    aux_phi_1 = sqrtf( - aux_sigma2sur4moinsa * psi_k + x_1 /expo ) + sigma / 2 * sqrtf(dt) * y;
                    x_1 = 1 / expo * aux_phi_1 * aux_phi_1 - aux_sigma2sur4moinsa * psi_k;
                }
                else {
                    u_tilde_1 = x_1 * expo_2 + aux_u_tilde_1;
                    u_tilde_2 = u_tilde_1 * u_tilde_1 + sigma * sigma * psi_k_2 * (a * psi_k_2 / 2.0f + x_1 * expo_2);
                    pi = 0.5f * (1 - sqrtf(1 - u_tilde_1 * u_tilde_1 / u_tilde_2) );
                    u_2 = curand_uniform(&localState);
                    if (u_2 < pi){
                        x_1 = u_tilde_1 / 2.0f / pi;
                    }
                    else {
                        x_1 = u_tilde_1 / 2.0f / (1.0f - pi);
                    };
                };
/////////////// Fin CIR_O2
                dx_1 += x_1;
                x_2 += (x_1 - 0.5*dx_1) * dt;
                x_4 += 0.5*x_3*dt;
                x_3 = x_3 * expf( aux_hw_1 + aux_hw_2 * dx_1 + aux_hw_3 * (x_1 - 0.5 * dx_1) );
                x_4 += 0.5*x_3*dt;
/////////////// HZ
                x_3 = x_3 * expf(sqrtf(x_1 * aux_hz_1)*n.x);
            }
        };
        if(type == 'e'){
            value_option = fmaxf(0, strike - x_3);
        }
        else {
            value_option = fmaxf(0, strike - x_4);
        }
        partial_sums[threadIdx.x] += value_option;
        partial_sums_squared[threadIdx.x] += value_option * value_option;
        // printf("%f\n", value_option);
    };

    // Synchronize the threads
    __syncthreads();

    // Sum per block
    if(threadIdx.x == 0){
        float sum = 0;
        float sum_squared = 0.0f;
        for (int i = 0; i < blockDim.x; ++i){
            sum += partial_sums[i];
            sum_squared += partial_sums_squared[i];
        };
        results_sum[blockIdx.x] += sum;
        results_sum_squared[blockIdx.x] += sum_squared;
    };

    // Saving the states in the global memory
    states[id] = localState;
};


long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
    return elapsed;
};


void wrapper_kernel_o2(float* output, float x_0, float cir_0, float r, float a, float k, float sigma, float rho, float expiry, float strike, unsigned int num_steps, char type){
    clock_t t1;
    clock_t t2;
    long elapsed;

    t1 = clock();
    float *h_results_sum, *d_results_sum;
    h_results_sum = (float*)malloc(num_blocks * sizeof(float));
    cudaMalloc(&d_results_sum, num_blocks * sizeof(float));
    for(int i = 0; i < num_blocks; ++i){
        h_results_sum[i] = 0;
    };
    float *h_results_sum_squared, *d_results_sum_squared;
    h_results_sum_squared = (float*)malloc(num_blocks * sizeof(float));
    cudaMalloc(&d_results_sum_squared, num_blocks * sizeof(float));
    for(int i = 0; i < num_blocks; ++i){
        h_results_sum_squared[i] = 0;
    };
    cudaMemcpy(d_results_sum, h_results_sum, num_blocks * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_results_sum_squared, h_results_sum_squared, num_blocks * sizeof(float), cudaMemcpyHostToDevice);

    curandState *d_states;
    cudaMalloc(&d_states, num_threads * num_blocks * sizeof(curandState));

    setup_states <<< num_blocks, num_threads >>> (d_states);
    cudaDeviceSynchronize();
    put_euro_heston_cir_o2 <<< num_blocks, num_threads >>> (d_states, d_results_sum, d_results_sum_squared, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, type);
    cudaDeviceSynchronize();

    cudaMemcpy(h_results_sum, d_results_sum, num_blocks * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results_sum_squared, d_results_sum_squared, num_blocks * sizeof(float), cudaMemcpyDeviceToHost);

    // Global sum
    float global_results_sum = 0;
    for(int i = 0; i < num_blocks; i++){
        global_results_sum += h_results_sum[i];
    };
    global_results_sum = global_results_sum * expf(-r * expiry);
    float global_results_sum_squared = 0;
    for(int i = 0; i < num_blocks; i++){
        global_results_sum_squared += h_results_sum_squared[i];
    };
    global_results_sum_squared = global_results_sum_squared * expf(- 2 * r * expiry);
    unsigned int num_simu = (num_blocks * num_threads * num_iterations);
    float empirical_expectency = global_results_sum / num_simu;
    float empirical_squared = global_results_sum_squared / num_simu;
    float empirical_variance = empirical_squared - empirical_expectency * empirical_expectency;
    float confidence_interval_low = empirical_expectency - 1.96 * sqrtf(empirical_variance / num_simu);
    float confidence_interval_high = empirical_expectency + 1.96 * sqrtf(empirical_variance / num_simu);

    t2 = clock();
    elapsed = timediff(t1, t2);

    output[0] = empirical_expectency;
    output[1] = confidence_interval_low;
    output[2] = elapsed;

    printf("We have computed a MC Call price of : %f\n", empirical_expectency);
    printf("Empirical variance : %f\n", empirical_variance);
    printf("Number of simulations: %d\n", num_simu);
    printf("Confidence interval : (%f , %f)\n", confidence_interval_low, confidence_interval_high);
    printf("Time elapsed: %ld ms\n", elapsed);

    return;
};


__global__ void put_heston_cir_o3(curandState* states, float* results_sum, float* results_sum_squared, float x_0, float cir_0, float r, float a, float k, float sigma, float rho, float t, float strike, unsigned int num_steps, char type)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    // Saving the state in the GPU memory to be more efficient
    curandState localState = states[id];
    // Shared memory for the Monte Carlo
    __shared__ float partial_sums[num_threads];
    __shared__ float partial_sums_squared[num_threads];

    // Some values that we will not have to compute twice / maybe we should give this in parameters.
    float x_1;
    float x_2;
    float x_3;
    float x_4;
    float dx_1;
    float dt = t / num_steps;
    float u_1; // u will be a uniform variable
    float u_2; // idem
    float u_3; // idem
    float2 n; // n will be a pair of normal variable
    float y;
    float epsilon;
    float u_tilde_1;
    float u_tilde_2;
    float u_tilde_3;
    float s;
    float p;
    float delta;
    float pi;
    float value_option;

    float sqrt3 = sqrtf(3.);
    float sqrt2 = sqrtf(2.);
    float s_3_m_s6 = sqrtf( 3 - sqrtf(6));
    float s_3_p_s6 = sqrtf( 3 + sqrtf(6));
    float aux_proba_y = (sqrtf(6.) - 2.)/ (4.0f*sqrtf(6.));

    float sigma_2 = sigma * sigma;
    float four_a_over_3 = 4*a / 3;
    float four_a = 4*a;
    float sigma_2_over_4_minus_a = sigma_2 / 4 - a;
    float sigma_2_over_4_minus_a_abs;
    if (sigma_2_over_4_minus_a > 0){
        sigma_2_over_4_minus_a_abs = sigma_2_over_4_minus_a;
    } else {
        sigma_2_over_4_minus_a_abs = -sigma_2_over_4_minus_a;
    };
    float aux_hz_1 = (1. - rho*rho) * dt;
    float aux_hw_1 = (r - rho * a / sigma) * dt;
    float aux_hw_2 = rho / sigma;
    float aux_hw_3 = (rho*k/sigma - 0.5) * dt;

    float psi_k;
    if (k == 0){
        psi_k = dt ;
    }
    else {
        psi_k = ( 1 - expf(-k * dt) ) / k;
    };
    float psi_minus_k;
    if (k == 0){
        psi_minus_k = dt;
    }
    else {
        psi_minus_k = ( 1 - expf(+k * dt) ) / (-k);
    };

    float k_3 = 0;
    float aux_k_3_1 = sqrtf(sigma_2 / 4 - a + sigma / sqrtf(2) * sqrtf(-sigma_2_over_4_minus_a)) + sigma / 2 * s_3_p_s6;
    float aux_k_3_2 = sigma / sqrt2 * sqrtf(-sigma_2_over_4_minus_a);
    float aux_k_3_3 = sqrtf(sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a)) + sigma / 2 * s_3_p_s6;
    float aux_k_3_4 = sigma_2 / 4 - a + aux_k_3_3 * aux_k_3_3;
    if (sigma_2 <= four_a_over_3) {
        k_3 = aux_k_3_2;
    } else if (sigma_2 <= four_a) {
        k_3 = aux_k_3_1 * aux_k_3_1;
    } else {
        k_3 = aux_k_3_4;
    };
    k_3 *= psi_minus_k;
    float expo = expf(k * dt / 2);
    float expo_2 = expf(-k*dt);
    float aux_u_tilde_1 = a*psi_k;
    float aux_u_tilde_3_1 = psi_k * (a + sigma_2 / 2) ;
    float aux_u_tilde_3_2 = 2*expf(-2*k*dt);

    // Main loop for the Monte Carlo
    for(unsigned int i = 0; i < num_iterations; ++i){
        // Initialization of the shared memory at the begining of the MC
        if(i == 0){
            partial_sums[threadIdx.x] = 0;
            partial_sums_squared[threadIdx.x] = 0;
        };
        x_1 = cir_0;  // Vol process
        x_2 = 0;  // Integration of the vol process
        x_3 = x_0;  // Stock process
        x_4 = 0;  // Integration of the stock process

        // CIR_O3 and Heston
        for (unsigned int k = 0; k < num_steps; ++k){
            // printf("Value of the Heston %f  and the CIR %f at the step %d. \n", x_3, x_1, k);
            u_1 = curand_uniform(&localState);
            n = curand_normal2(&localState); // It is not optimal, we simulate two uniform for 1 normal
            if (u_1 < 0.5) {
                // HZ
                x_3 = x_3 * expf(sqrtf(x_1 * aux_hz_1)*n.x);
                // HW
                dx_1 = - x_1;
/////////////// CIR_O3
                if(x_1 >= k_3){
                    u_1 = curand_uniform(&localState);
                    u_2 = curand_uniform(&localState);
                    u_3 = curand_uniform(&localState);
                    // Computing y
                    if (u_1 < aux_proba_y) {
                        y = - s_3_p_s6;
                    } else if (u_1 < 2. * aux_proba_y){
                        y = + s_3_p_s6;
                    } else if (u_1 < 0.5 + aux_proba_y){
                        y = - s_3_m_s6;
                    } else {
                        y = + s_3_m_s6;
                    };
                    // Computing epsilon
                    if (u_2 < 1./2.){
                        epsilon = -1.;
                    } else {
                        epsilon = +1.;
                    };
                    // zeta
                    if (u_3 < 1./3.) {
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                        } else{
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                        };
                    } else if (u_3 < 2./3.){
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2.0f) * fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                        } else{
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 = fmaxf(0.0f,  sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2.0f) * fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                        };
                    } else {
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 = fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2.f) * fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                        } else{
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 = fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2.0f) * fmaxf(0.0f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                        };
                    };
                    x_1 *= expo_2;
                } else {
                    u_1 = curand_uniform(&localState);
                    u_tilde_1 = x_1 * expo_2 + aux_u_tilde_1;
                    u_tilde_2 = u_tilde_1 * u_tilde_1 + sigma * sigma * psi_k * (a * psi_k  / 2.0f + x_1 * expo_2);
                    u_tilde_3 = u_tilde_1 * u_tilde_2 + sigma_2 * psi_k * (x_1*x_1 * aux_u_tilde_3_2 + aux_u_tilde_3_1* (3*x_1 * expo_2 + a * psi_k));
                    s = (u_tilde_3 - u_tilde_1 * u_tilde_2) / (u_tilde_2 - u_tilde_1 * u_tilde_1);
                    p = (u_tilde_1 * u_tilde_3 - u_tilde_2*u_tilde_2) / (u_tilde_2 - u_tilde_1*u_tilde_1);
                    delta = sqrtf(s*s - 4.*p);
                    pi = (u_tilde_1 - (s-delta) / 2.f) / delta;
                    if (u_1 < pi) {
                        x_1 = (s + delta) / 2.f;
                    } else {
                        x_1 = (s - delta) / 2.f;
                    };
                };
/////////////// Fin CIR_O3
                dx_1 += x_1;
                x_2 += (x_1 - 0.5*dx_1) * dt;
                x_4 += 0.5*x_3*dt;
                x_3 = x_3 * expf( aux_hw_1 + aux_hw_2 * dx_1 + aux_hw_3 * (x_1 - 0.5 * dx_1) );
                x_4 += 0.5*x_3*dt;
            }
            else {
                // HW
                dx_1 = - x_1;
/////////////// CIR_O3
                if(x_1 >= k_3){
                    u_1 = curand_uniform(&localState);
                    u_2 = curand_uniform(&localState);
                    u_3 = curand_uniform(&localState);
                    // Computing y
                    if (u_1 < aux_proba_y) {
                        y = - s_3_p_s6;
                    } else if (u_1 < 2. * aux_proba_y){
                        y = + s_3_p_s6;
                    } else if (u_1 < 0.5 + aux_proba_y){
                        y = - s_3_m_s6;
                    } else {
                        y = + s_3_m_s6;
                    };
                    // Computing epsilon
                    if (u_2 < 1./2.){
                        epsilon = -1.;
                    } else {
                        epsilon = +1.;
                    };
                    // zeta
                    if (u_3 < 1./3.) {
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 = fmaxf(0.f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2.f) * fmaxf(0.f, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                        } else{
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                        };
                    } else if (u_3 < 2./3.){
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                        } else{
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                        };
                    } else {
                        if (sigma_2_over_4_minus_a <= 0) {
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                        } else{
                            x_1 += sigma / sqrt2 * sqrtf(sigma_2_over_4_minus_a_abs) * epsilon * psi_minus_k;  // Xt
                            x_1 += -sigma_2_over_4_minus_a * psi_minus_k;  // X0
                            x_1 = fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2) * fmaxf(0, sqrtf(x_1) + sigma*sqrtf(psi_minus_k)*y/2); // X1
                        };
                    };
                    x_1 *= expo_2;
                } else {
                    u_1 = curand_uniform(&localState);
                    u_tilde_1 = x_1 * expo_2 + aux_u_tilde_1;
                    u_tilde_2 = u_tilde_1 * u_tilde_1 + sigma * sigma * psi_k  * (a * psi_k  / 2.0f + x_1 * expo_2);
                    u_tilde_3 = u_tilde_1 * u_tilde_2 + sigma_2 * psi_k * (x_1*x_1 * aux_u_tilde_3_2 + aux_u_tilde_3_1* (3*x_1 * expo_2 + a * psi_k));
                    s = (u_tilde_3 - u_tilde_1 * u_tilde_2) / (u_tilde_2 - u_tilde_1 * u_tilde_1);
                    p = (u_tilde_1 * u_tilde_3 - u_tilde_2*u_tilde_2) / (u_tilde_2 - u_tilde_1*u_tilde_1);
                    delta = sqrtf(s*s - 4.*p);
                    pi = (u_tilde_1 - (s-delta) / 2.) / delta;
                    if (u_1 < pi) {
                        x_1 = (s + delta) / 2;
                    } else {
                        x_1 = (s - delta) / 2;
                    };
                };
/////////////// Fin CIR_O3
                dx_1 += x_1;
                x_2 += (x_1 - 0.5*dx_1) * dt;
                x_4 += 0.5*x_3*dt;
                x_3 = x_3 * expf( aux_hw_1 + aux_hw_2 * dx_1 + aux_hw_3 * (x_1 - 0.5 * dx_1) );
                x_4 += 0.5*x_3*dt;
/////////////// HZ
                x_3 = x_3 * expf(sqrtf(x_1 * aux_hz_1)*n.x);
            }
        };
        if(type == 'e'){
            value_option = fmaxf(0, strike - x_3);
        }
        else {
            value_option = fmaxf(0, strike - x_4);
        };

        partial_sums[threadIdx.x] += value_option;
        partial_sums_squared[threadIdx.x] += value_option * value_option;
        // printf("%f\n", value_option);
    };

    // Synchronize the threads
    __syncthreads();

    // Sum per block
    if(threadIdx.x == 0){
        float sum = 0;
        float sum_squared = 0.0f;
        for (int i = 0; i < blockDim.x; ++i){
            sum += partial_sums[i];
            sum_squared += partial_sums_squared[i];
        };
        results_sum[blockIdx.x] += sum;
        results_sum_squared[blockIdx.x] += sum_squared;
    };

    // Saving the states in the global memory
    states[id] = localState;
};

void wrapper_kernel_o3(float* output, float x_0, float cir_0, float r, float a, float k, float sigma, float rho, float expiry, float strike, unsigned int num_steps, char type){
    clock_t t1;
    clock_t t2;
    long elapsed;

    t1 = clock();
    float *h_results_sum, *d_results_sum;
    h_results_sum = (float*)malloc(num_blocks * sizeof(float));
    cudaMalloc(&d_results_sum, num_blocks * sizeof(float));
    for(int i = 0; i < num_blocks; ++i){
        h_results_sum[i] = 0;
    };
    float *h_results_sum_squared, *d_results_sum_squared;
    h_results_sum_squared = (float*)malloc(num_blocks * sizeof(float));
    cudaMalloc(&d_results_sum_squared, num_blocks * sizeof(float));
    for(int i = 0; i < num_blocks; ++i){
        h_results_sum_squared[i] = 0;
    };
    cudaMemcpy(d_results_sum, h_results_sum, num_blocks * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_results_sum_squared, h_results_sum_squared, num_blocks * sizeof(float), cudaMemcpyHostToDevice);

    curandState *d_states;
    cudaMalloc(&d_states, num_threads * num_blocks * sizeof(curandState));

    setup_states <<< num_blocks, num_threads >>> (d_states);
    cudaDeviceSynchronize();
    put_heston_cir_o3 <<< num_blocks, num_threads >>> (d_states, d_results_sum, d_results_sum_squared, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, type);
    cudaDeviceSynchronize();

    cudaMemcpy(h_results_sum, d_results_sum, num_blocks * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_results_sum_squared, d_results_sum_squared, num_blocks * sizeof(float), cudaMemcpyDeviceToHost);

    // Global sum
    float global_results_sum = 0;
    for(int i = 0; i < num_blocks; i++){
        global_results_sum += h_results_sum[i];
    };
    global_results_sum = global_results_sum * expf(-r * expiry);
    float global_results_sum_squared = 0;
    for(int i = 0; i < num_blocks; i++){
        global_results_sum_squared += h_results_sum_squared[i];
    };
    global_results_sum_squared = global_results_sum_squared * expf(- 2 * r * expiry);
    unsigned int num_simu = (num_blocks * num_threads * num_iterations);
    float empirical_expectency = global_results_sum / num_simu;
    float empirical_squared = global_results_sum_squared / num_simu;
    float empirical_variance = empirical_squared - empirical_expectency * empirical_expectency;
    float confidence_interval_low = empirical_expectency - 1.96 * sqrtf(empirical_variance / num_simu);
    float confidence_interval_high = empirical_expectency + 1.96 * sqrtf(empirical_variance / num_simu);

    t2 = clock();
    elapsed = timediff(t1, t2);

    output[0] = empirical_expectency;
    output[1] = confidence_interval_low;
    output[2] = elapsed;

    printf("We have computed a MC Call price of : %f\n", empirical_expectency);
    printf("Empirical variance : %f\n", empirical_variance);
    printf("Number of simulations: %d\n", num_simu);
    printf("Confidence interval : (%f , %f)\n", confidence_interval_low, confidence_interval_high);
    printf("Time elapsed: %ld ms\n", elapsed);

    return;
};

void cuda_plot_graph_performance(double expiry, double strike, double cir_0, double x_0, double a, double k , double sigma, double rho,double r, char type, int* num_steps_array, int num_points, unsigned int cap, double precision, double exact_value){

    float output_2[3];
    float output_3[3];

    FILE *f = fopen("plot.dat", "w");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    };
    unsigned int num_steps;
    for (unsigned int i = 0; i < num_points; ++i){
        num_steps = num_steps_array[i];
        wrapper_kernel_o2(output_2, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, type);
        wrapper_kernel_o3(output_3, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, type);
        fprintf(f, "%f\t%f\t%f\t%f\t%f\n", 1./num_steps, output_2[0], output_2[1], output_3[0], output_3[1]);
    };
    printf("We have computed the trajectories.\n");
    fclose(f);
    printf("We have written the data.");
    *f = fopen("gnu", "w");
    fprintf(f,"set nokey\n") ;
    fprintf(f, "set xlabel \"Inverse of number of steps\"\n");
    fprintf(f, "plot ");
    fprintf(f, "\"plot.dat\" using 1:2:3 with yerrorlines, \\\n");
    fprintf(f, "\"plot.dat\" using 1:4:5 with yerrorlines, \\\n");
    fprintf(f, "with lines lt 3");
    fclose(f);

    printf("We have written the gnuplot file.\n");

    return;
};

int main(void){
    float k = 0.5f;
    float a = 0.02f;
    float sigma = 0.4f;
    float x_0 = 100.0f;
    float cir_0 = 0.04;
    float rho = - 0.5f;
    float r = 0.02f;
    float strike = 100.0f;
    float expiry = 1.0f;
    unsigned int num_steps = 100;

    float output[3];
    wrapper_kernel_o2(output, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, 'e');
    wrapper_kernel_o3(output, x_0, cir_0, r, a, k, sigma, rho, expiry, strike, num_steps, 'e');

    printf("%f\t%f\t%f", output[0], output[1], output[2]);
    return 0;
}
