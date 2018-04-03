#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>

#include <time.h>

const int num_blocks = 1024;
const int num_threads = 256;
const int num_iterations = 10;
const int num_steps = 100;

__global__ void setup_states(curandState* states){
    int id = threadIdx.x + num_threads * blockIdx.x;
    // Initialisation states
    curand_init(0, id, 0, &states[id]);
}

__global__ void put_euro_heston_cir_o2(curandState* states, float* results_sum, float* results_sum_squared, float x_0, float cir_0, float r, float a, float k, float sigma, float rho, float t, float strike)
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
        value_option = fmaxf(0, strike - x_3);
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


int main(void){

    clock_t t1;
    clock_t t2;
    long elapsed;

    float k = 0.5f;
    float a = 0.02f;
    float sigma = 0.4f;
    float x_0 = 100.0f;
    float cir_0 = 0.04;
    float rho = - 0.5f;
    float r = 0.02f;
    float strike = 100.0f;
    float expiry = 1.0f;


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
    put_euro_heston_cir_o2 <<< num_blocks, num_threads >>> (d_states, d_results_sum, d_results_sum_squared, x_0, cir_0, r, a, k, sigma, rho, expiry, strike);
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

    printf("We have computed a MC Call price of : %f\n", empirical_expectency);
    printf("Empirical variance : %f\n", empirical_variance);
    printf("Number of simulations: %d\n", num_simu);
    printf("Confidence interval : (%f , %f)\n", confidence_interval_low, confidence_interval_high);
    elapsed = timediff(t1, t2);
    printf("Time elapsed: %ld ms\n", elapsed);

    return 0;
}
