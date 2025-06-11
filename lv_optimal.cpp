//  file: lv_optimal.cpp
//
//  C++ Program to demonstrate optimize the parameters of the lv model
//
//  Programmer:  Jake Kamen jkamen23@osu.edu
//
//  Revision history:
//      03/24/25  original version, based on discussion searches of 
//                the interwebs
//      04/16/25  comments added for easier reading and understanding
//
//  Notes:  
//   * Compile and link with:
//       make -f make_lv_optimal
//   * A much more in-depth description of this file and the process are
//     included in the Overview.cpp file
//   * All sections until the end of the 4th order runge kutta are essentially
//     the same as the lotka_volterra.cpp file, with additional #includes
//   * The fitness function isn't great and a more refined version should be 
//     built if work on this project continues
// 
//*********************************************************************//
// Retrive necessary inclusions
#include <iostream>  // for input and output operations
#include <vector>    // for storing dynamic arrays of parameters and results
#include <tuple>     // for handling multiple returned values (in rk4)
#include <iomanip>   // fomatting outputs
#include <fstream>   // for file operation
#include <cmath>     // for math functions
#include <limits>    // for extreme values

using namespace std;

// fmt will be used in output formatting
#define fmt setw(20) << setprecision(10) << scientific

// Parameters for Lotka-Volterra
double alpha, beta, delta, eta, h;

// Lotka-Volterra equations
// For prey
double dprey(double prey, double predator, double alpha_val, double beta_val) {
    return (prey * alpha_val - prey * predator * beta_val);
    //       prey growth       prey reduction from predation
}
// For predators
double dpredator(double prey, double predator, double delta_val, double eta_val) {
    return (predator * delta_val * prey - eta_val * predator);
    //      predator growth from eating    predator death from starvation
}

// A 4th-order Runge-Kutta method
auto rk4(double prey, double predator, double alpha_val, double beta_val, double delta_val, double eta_val) {
    double kx1 = h * dprey(prey, predator, alpha_val, beta_val);
    double ky1 = h * dpredator(prey, predator, delta_val, eta_val);  
    double kx2 = h * dprey(prey + 0.5 * kx1, predator + 0.5 * ky1, alpha_val, beta_val);
    double ky2 = h * dpredator(prey + 0.5 * kx1, predator + 0.5 * ky1, delta_val, eta_val);
    double kx3 = h * dprey(prey + 0.5 * kx2, predator + 0.5 * ky2, alpha_val, beta_val);
    double ky3 = h * dpredator(prey + 0.5 * kx2, predator + 0.5 * ky2, delta_val, eta_val);
    double kx4 = h * dprey(prey + kx3, predator + ky3, alpha_val, beta_val);
    double ky4 = h * dpredator(prey + kx3, predator + ky3, delta_val, eta_val);
    
    return make_tuple(prey + (kx1 + 2 * (kx2 + kx3) + kx4) / 6,
                      predator + (ky1 + 2 * (ky2 + ky3) + ky4) / 6);
}

// Fitness function based to check for oscillatory behavior
double fitness_function(double alpha_val, double beta_val, double delta_val, double eta_val) {

    // Parameters set outside of the 4 in the equation, such as initial 
    // starting population sizes. Prey and predator minimums and maximums
    // are also included for the changing parameters
    double prey = 100, predator = 10, time = 0.0, endtime = 50.0;
    double prey_max = 0.0, prey_min = 1e6, predator_max = 0.0, predator_min = 1e6;
    
    // Step size
    h = 0.1; 

    // Simulate over time
    while (time <= endtime) {        
        tie(prey, predator) = rk4(prey, predator, alpha_val, beta_val, delta_val, eta_val);
        
       
        // Track the maximum and minimum values for prey and predator
        prey_max = max(prey_max, prey);
        prey_min = min(prey_min, prey);
        predator_max = max(predator_max, predator);
        predator_min = min(predator_min, predator);

        time += h;
    }

    // Calculate the amplitudes for prey and predator populations
    double prey_amplitude = prey_max - prey_min;
    double predator_amplitude = predator_max - predator_min;

    // Return a combined measure of both prey and predator amplitudes
    // This will help us to tell if this oscillatory motion exists
    return (prey_amplitude + predator_amplitude);
}

// Now that all the functions have been defined, we can begin to run it across a range of parameters

int main() {
    double best_fitness = std::numeric_limits<double>::infinity();
    vector<double> best_params;

    // Define ranges and step sizes for each parameter
    double alpha_min = 0.1, alpha_max = 1.0, alpha_step = 0.1;
    double beta_min = 0.001, beta_max = 0.1, beta_step = 0.001;
    double delta_min = 0.001, delta_max = 0.1, delta_step = 0.001;
    double eta_min = 0.01, eta_max = 1.0, eta_step = 0.01;

    // Perform search over the parameter space
    for (double alpha_val = alpha_min; alpha_val <= alpha_max; alpha_val += alpha_step) {
    cout << "We are on step:" << alpha_val << endl;
        for (double beta_val = beta_min; beta_val <= beta_max; beta_val += beta_step) {
            for (double delta_val = delta_min; delta_val <= delta_max; delta_val += delta_step) {
                for (double eta_val = eta_min; eta_val <= eta_max; eta_val += eta_step) {
                    // Evaluate fitness for current parameter combination
                    double fitness = fitness_function(alpha_val, beta_val, delta_val, eta_val);

                    // Track the best combination (lowest fitness value)
                    if (fitness < best_fitness) {
                        best_fitness = fitness;
                        best_params = {alpha_val, beta_val, delta_val, eta_val};
                    }
                }
            }
        }
    }

    // Output the optimal parameters
    cout << "Optimal Parameters: \n";
    cout << "Alpha: " << best_params[0] << " Beta: " << best_params[1] 
         << " Delta: " << best_params[2] << " Eta: " << best_params[3] << endl;

    return 0;
}
