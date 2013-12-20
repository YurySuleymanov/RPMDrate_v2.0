#!/usr/bin/env python
# encoding: utf-8
  
from PES import get_potential, initialize_potential 
  
################################################################################ 
  
label =  'Cl + O3 -> O2 + ClO' 
  
reactants( 
    atoms = ['O', 'O', 'O', 'Cl'], 
    reactant1Atoms = [1,2,3], 
    reactant2Atoms = [4], 
    Rinf = (8,"angstrom"), 
) 
  
transitionState( 
    geometry = ( 
         [[-0.79756700, 3.40691242, 0.61271692], 
         [-0.24667794, 2.64101070, -0.30298481], 
         [1.04797818, 2.66283161, -0.30298481], 
         [1.89021513, 1.18915913, 1.46151001]], 
        "angstrom"), 
    formingBonds = [(3,4)], 
    breakingBonds = [(2,3)], 
) 
equivalentTransitionState( 
    formingBonds=[(1,4)],  
    breakingBonds=[(2,1)], 
) 
  
thermostat('Andersen') 
  
################################################################################# 
  
xi_list = numpy.arange(-0.05, 1.11, 0.01) 
#xi_list = [1.0] 
   
generateUmbrellaConfigurations( 
    dt = (0.0001,"ps"), 
    evolutionTime = (5,"ps"), 
    xi_list = xi_list, 
    kforce = 0.1 * T, 
)  
   
#xi_list = numpy.arange(-0.05, 1.11, 0.01) 
#xi_list = [-0.05]
 
windows = [] 
for xi in xi_list: 
    window = Window(xi=xi, kforce=0.1*T, trajectories=118, equilibrationTime=(20,'ps'), evolutionTime=(100,'ps'))  
    windows.append(window) 
 
conductUmbrellaSampling( 
    dt = (0.0001,"ps"), 
    windows = windows, 
) 
 
computePotentialOfMeanForce(windows=windows, xi_min=-0.04, xi_max=1.10, bins=5000) 
 
computeRecrossingFactor( 
    dt = (0.0001,"ps"), 
    equilibrationTime = (20,"ps"), 
    childTrajectories = 23900, 
    childSamplingTime = (2,"ps"), 
    childrenPerSampling = 100, 
    childEvolutionTime = (0.30,"ps"), 
) 
   
computeRateCoefficient() 
