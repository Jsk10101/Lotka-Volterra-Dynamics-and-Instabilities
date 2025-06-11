//  file: lotka_volterra.cpp
//
//  C++ Program to demonstrate the use of the basic lotka volterra model
//
//  Programmer:  Jake Kamen jkamen23@osu.edu
//
//  Revision history:
//      03/24/25  original version, based on discussion in Youtube videos
//                and equations given within databases on the internet
//      04/16/25  Revised with real world scenarios and necessary commentation
//
//  Notes:  
//   * Compile and link with:
//       make -f make_lotka_volterra
//   * For a much more in-depth description of the development of this .cpp
//     file as well as the resulting makefile, .x file, .dat file, and 
//     .plt file, look at the Overview.cpp document, specifically 
//     "Section 2 - lotka_volterra.cpp"
//   * If you a wondering about parameter choices, look at "Section 3 - Real
//     World Scenarios" and "Section 4 - lv_optimal.cpp" of Overview.cpp
//   * Definetly check out "Section 3 - Real World Scenarios" to see how
//     this model can be applied to valuable environmental situations
//
//*********************************************************************//
// Retrive necessary inclusions
#include<tuple>
#include<iomanip>
#include <iostream>
#include<fstream>
using namespace std;

// fmt will be used in output formatting
#define fmt  setw(20) << setprecision(10) << scientific 

// define parameters that will be used as well as step size
double alpha, beta, delta, eta, h;

// differential equation for the prey population
double dprey(double prey,double predator){
	return(prey*alpha - prey*predator*beta);
}

// differential equation for the predator population
double dpredator(double prey,double predator){
	return(predator*delta*prey - eta*predator);
}

// A 4th-order Runge-Kutta method
auto rk4(double prey,double predator){
	//K1
	double kx1 = h*dprey(prey,predator);
	double ky1 = h*dpredator(prey,predator);	
	//k2
	double kx2 = h*dprey(prey+.5*kx1,predator+.5*ky1);
	double ky2 = h*dpredator(prey+.5*kx1,predator+.5*ky1);
	//K3
	double kx3 = h*dprey(prey+.5*kx2,predator+.5*ky2);
	double ky3 = h*dpredator(prey+.5*kx2,predator+.5*ky2);
	//K4
	double kx4 = h*dprey(prey+kx3,predator+ky3);
	double ky4 = h*dpredator(prey+kx3,predator+ky3);
	
	return make_tuple(prey+(kx1+2*(kx2+kx3)+kx4)/6,predator+(ky1+2*(ky2+ky3)+ky4)/6);
}

// Now we can run the program
int main(){

	//System parameters
	double prey=110, predator=12,time=0.0,endtime=10.0;		
	
	//step size
	h     = 1e-4;
	
	// These parameters below are the parameters derived from the final
	// run of "Section 4 - lv_optimal.cpp" in Overview.cpp. To learn more
	// about why they work, look at that file
	alpha = 0.1; 
	beta  = 0.01; 
	delta = 0.001;
        eta   =	0.1;
	
	/*
	Alpha : Natural Birth rate of prey in absence of any predation
	Beta  : Death rate of preys due to predation (hunting efficiency)
        eta   : Death rate of predators in the absence of prey
	Delta : Rate that predators increase due to consuming prey
	h     : step-size    
	*/
	
	/*
	Differential equation used for numerical simulations.
	prey'	 	= alpha*prey - beta*prey*predator
	predator' 	= delta*prey*predator*beta - eta*predator
	*/
	
	// output file for later
	ofstream output("lotka_volterra.dat");
	
	output <<"#"<< setw(19) <<"Time" << setw(20) <<"Prey" << setw(20) <<"Predator" <<endl;
	while(time<=endtime){		
		tie(prey,predator)=rk4(prey,predator);
		output << fmt << time << fmt << prey << fmt << predator << endl;
		time+=h;
	}
	output.close();	
	
	// Write this to make sure the code runs and data has been saved
	cout << "The data has been saved to lotka_volterra.dat" << endl;
	
}
