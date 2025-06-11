# Lotka-Volterra-Dynamics-and-Instabilities
This document is simply an overview of this entire project. Since no other document is sumbitted with this project, it will serve as a tell all document for each section, describing code, data collection, and implementation. If anything is not included in the file itself or is discussed sparsely, it will be included here in much greater detail.


Section 1 - The Lotka Volterra Model
My project focuses on the Lotka-Volterra Model, a pair of first order nonlinear differential equations that describe the interactions between predators and prey. The idea for the project is to start fairly simply, making sure this basic system works in code, and then apply more complex scenarios and stipulations to see how it changes outcomes. To start, I will try to implement the basic model into code done through the lotka_volterra.cpp file. Online searching has gotten us the coupled equations and their parameters, which will be specified in the .cpp documents.


Section 2 - lotka_volterra.cpp
For the differential equation for prey, it comes in the form:
dprey = (prey*alpha - prey*predator*beta),
where alpha describes the natural birth rate of prey in absence of any predation and beta describes the death rate of prey due to predation, a variable called hunting efficiency. This dprey value simply describes the rate of change in the prey population, acting as the stand in for the differential of prey population over time
For the differential equation for predators, it is in the form:
dpredators = (predator*delta*prey - eta*predator),
where delta describes the rate that predators increase due to consuming prey and eta describes the death rate of predators in the absence of prey.

To tie the two equations together and integrate, we use the 4th order Runge-Kutta method, which allows us to discretize the approximate solutions for the nonlinear equations. Once set up, we can deal with our 'main' region of the code, setting our initial prey and predator population sizes as well as setting values for the 4 parameters described above. The parameters are set internally, not set once the .x file is run, to help me keep track of their changing values. Since much of this project is slighty changing these, its easier to remember where I'm at if I keep the values saved in the .cpp file. For the initial case, we will use parameters set for a well-known case, where alpha = 1.1, beta = 0.4, eta = 0.4, and delta = 0.1. This is out initial validation set, making sure that the lotka volterra model works generally. To save the output values for later, an ofstream output is written to the file "lotka_volterra.dat". This file saves the prey and predator populations at each iterative step to allow us to plot populations as time changes later on.

Now that the .cpp file has been written, we need a make file to run it. This file, called "make_lotka_volterra" allows us to run the "lotka_volterra.cpp" file, which creates our .x file. To write the make file, I simply used previous makefiles that have worked in this class, not changing anything but the BASE and SRCS lines to match my file names.

Once the makefile is run, the lotka_volterra.x file can be run, providing us with our lotka_volterra.dat file that we will need to plotting. Taking that .dat file, we can use gnuplot to see our results. It's plotted using a .plt file, which makes labeling axis and columns for the data much easier. We can simply load in the .plt file in gnuplot and it gives us an image with good bounds, a title, and all the other things necessary for a complete plot. The file saved from this run is called "lv_a1.1_b0.4_c0.4_d0.1.ps" so that I you can easily know the value of the parameters input for the simulation.



Section 3 - Real World Scenarios
In this section, we will use the lotka_volterra.cpp and known parameters in the environment to simulate real populations. A wide range of real world scenarios have been gotten from the internet to see how a range of parameters affect the outcome of the lotka-volterra model. Below are a number of these relationships, the parameters that were gotten, and descriptions on if they worked out. If they did, the plots has been saved, the names of which you can find in each sub-section

Rabbits + Foxes (Australia):
Parameters - alpha = 1.0, beta = 0.1, delta = 0.075, eta = 1.5
            Rabbit Population = 40, Fox Population = 9
Produces a very quickly changing yet oscillatory relationship. This file has been plotted and saved as "lv_foxes_n_rabbits.ps" This relationship is particularly interesting because of how fast the populations can change. You can see that the eta value is extremely high compared to our first case done in section 2. In that case, about 4 cycles of oscillations occur in the timeframe of 50. In this rabbit and fox case, over 9 cycles happen. This timescale is an extremely noteworthy feature of this prey-predator relationship. Make note of this to compare with future animal relationship cycles.

Moose + Wolves (Isle Royale, MI):
Parameters - alpha = 0.1, beta = 0.005, delta = 0.0002, eta = 0.015
            Moose Population = 100, Wolf Population = 18
To get moose and wolf initial populations so that it presents oscillatory motion, we manually solve the differential equations for x and y, being the prey and predator populations, the resulting equilibrium populations come out to be x = eta/delta and y = alpha/beta. So by using our given parameters for alpha, beta, delta, and eta initially, we can solve to see the equilibrium population. We get moose=100 and wolves=20. We don't want exactly equilibrium, however, as we will not get oscillatory motion. Instead, we start the wolf population slightly below equilibrium, in this case to wolves=18. The resulting simulation is great, producing the type of motion we want. It can be seen in the "lv_wolves_n_moose.ps" file. It is very different visually to the fox and rabbit population. While those populations change very drastically up and down across the 50 time frame, the moose and wolf populations barely change over a 200 step time frame. Another key point is that, while the fox population drastically increase due to predation, the wolf population barely increases in the long run. Neither does the moose population, with both staying fairly close to their initial population values. This is a very different relationship to even our initial tested case, by which we know has well simulated oscillatory motion.

Phytoplankton + Zooplankton:
Parameters - alpha = 2.0, beta = 0.1, delta = 0.075, eta = 0.5
             Phytoplankton Population = 100, Zooplankton Population = 5
This scenario is definetly the oddest of the bunch. It helps to show, however, that this relationship spans across many species in many places. There are a few key differences. Firstly, the timescale for this plot are not in years like the other animal plots but are in days. These species reproduce unbelievably compared to other animals and this explains the drastic difference in timescale. Additionally, the units of population are not as singular like they would be for the other animals desribed. You obviously cannot have 1/2 of a moose in the real world. But the units for the plankton's are actually in biomass terms, so in g/m^3 terms. This means that a more continuous population change, as shown in the plots, actually have a more physical meaning for this scenario.

Throughout this section, there are a few pieces of complexity overlooked that would further refine a real world model. Mainly, a population in a specific environment often has a population cap. An environment can only hold so much of a species and so, when the population gets too high and there's not enough food to get around, the populations often begin to drop rapidly. There are other simplifications that have been added to the models. For example, these simulations assume closed systems, meaning there are not other predators that eat the prey, no other prey for the predators to eat, or other external factors such as weather. Our lotka-volterra model gives us a good idea of the relationships that populations have, but they still oversimplify them.



Section 4 - lv_optimal.cpp
lv_optimal.cpp tries to solve a new issue related to the lotka-volterra model: How do we find sets of parameters that actually work? It's great to find scenarios online that work in the real world, but if I had no prior knowledge to what works, what would I do to find a set that produces oscillatory motion? These are answered in this section of the project. Firstly, the file lv_optimal.cpp is extremely similar to the lotka_volterra.cpp file. It implements the same pair of differntial equations and again uses the 4th order runge kutta method to numerically solve. The new piece that is included here is the fitness function, a function defined in the code that checks the differences in the predator and prey populations across time. By seeing how these population values fluctuate, especially relative to one another, we can see if a model produced by some set of parameters gives oscillatory behavior or if it runs off one way or the other. Here, the initial population size must be set outside of these parameter checking loops. These components could be added on for additional layers of looping, getting a parameter space over all 6 input parameters for the model, but this would make running the code and getting good results unbelievably harder. As such, the parameters for initial population are set first. Then, a range from minimum to maximum alpha, beta, delta, and eta are set within the code, alongside a stepsize for each as it rotates through the for loop. If the resulting set of parameters are a best fit for the oscillatory motion, they are saved. At the end, whichever parameter space produces the best oscillations will be spit out into the bash cell.

Similar again to the lotka_volterra.cpp file, the lv_optimal.cpp file has its own make file called mave_lv_optimal. Again, not much in this file is changed from what we have done within the class except for the BASE and SRCS lines, where this file name is input. One run, a .x file is spit out. Once this file is run (which usually takes a while with a large parameter space), we get the output of the optimal set of parameters to produce oscillatory behavior.

Below is results for lv_optimal runs. They are broken up into what the initial population parameters are and then the resulting equational parameters, with the parenthesis after them representing the range of values from min to max that each variable was allowed to sweep through:

4a.
For populations - Prey = 20, Predator = 1
Result 1: alpha = 0.1 (0.1 - 2.0)
          beta = 0.1 (0.01 - 1.0)
          delta = 0.01 (0.01 - 0.3)
          eta = 0.1 (0.1 - 1.0)
          
Result 2: alpha = 0.01 (0.01 - 0.2)
          beta = 0.01 (0.001 - 0.05)
          delta = 0.001 (0.001 - 0.05)
          eta = 0.01 (0.01 - 0.1)

Result 3: alpha = 0.1 (0.1 - 1.0)
          beta = 0.1 (0.001 - 0.1)
          delta = 0.1(0.001 - 0.1)
          eta = 1.0 (0.01 - 1.0)
          
Result 3 has the widest range so most likely gives the best result. When plotting, we get a suprisingly nice oscillatory model, shown in file "lv_optimized_prey20_pred1.ps". One interesting feature that occurs is that, as described in the moose and wolf section, that the equilibrium occurs when       x = eta/delta and y = alpha/beta. For this result, the closest we can get equilibrium given the intital populations are the results that have been spit out. For result 3, alpha/beta = 0.1/0.1 = 1, which is the initial predator population. Using the same process, eta/delta = 1.0/0.1 = 10, which is the closest the prey population can get to its initial population 20 given the sweep of initial parameters. Let's try this again with a larger set of population, now accounting for these equilibrium relationships

4b.
For populations - Prey = 100, Predator = 10
Results: alpha = 0.1 (0.1 - 1.0)
          beta = 0.01(0.001 - 0.1)
          delta = 0.001 (0.001 - 0.1)
          eta = 0.1 (0.01 - 1.0)
          
The ranges of each parameter for this case has now been set so that it can reach a predator population a number of ways, with alpha/beta = 1.0/0.1 = 0.1/0.01 = 10. The same can be said for prey population, with eta/delta = 1.0/0.01 = 0.1/0.001 = 100. Giving the parameter space many ways to achieve the equalibrium population let it find one with this needed oscillary behavior. The biggest issue here is that it found a position in equilibrium. In 4a, preventing it from reaching a place of pure equilibrium meant it produced oscillatory behavior. In order to make sure this result plots and produces some form of oscillatory behavior, instead of pure equilibrium, we must tweak the intial populations slightly. So with alpha = 0.1, beta = 0.01, delta - 0.001, and eta = 0.1, we will tweak the initial population to prey = 110 and predator = 12. This model is then plotted, producing the plot "lv_optimized_prey100_pred10". It is important to remember that the pure equilibrium found and the population in this lv_optimal run is not the exact same population that is used to run the lotka_volterra.cpp file. 



Section 5 - Parameter Mapping
While this section could be much longer and have entire code chunks checking the results, I will have to leave it purely to discussion for now. Regardless, understanding the relationships that each parameter has throughout this notebook is extremely important. The relationship works both ways. Having equilibrium starting populations for the set of parameters in the differential equations means that if we know the parameters of the equations, we can set the the starting population so that it oscillates. If we instead know the starting population, we can solve for our equilibrium parameters to create oscillation. Alpha and beta are thus related to predator population and delta and eta to prey population. If alpha and beta are such that their result is close to the predator starting population, given the prey population will be in equilibrium, then the entire system will be in equilibrium. If delta and eta are such that their result is close to the prey starting population, given that the predator population will be in equilibrium, then the system will again be in euqilibrium. It's important that both of these scenarios are true at the same time, that both ratios are similar to their starting populations. They can be multiplicitavely different form the starting populations, as long as both the ratios are multiplicatively different from the starting populations by the same amount. So multiplying the both ratios by 2 to get a new starting population for both prey and predator will also produce an oscillatory scenario. In the future, acutal code to test this theory (or back it up) can be done.
