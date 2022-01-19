#!/usr/bin/python3
#
# A small program for investigating how the pheromone concentration
# of the original ACO algorithm
#
# Eero Holmstrom, 2019
#

import sys
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import parameters

# Usage
if len(sys.argv) != 6:
    print("Usage: pheromone_tester.py [tau_0] [rho] [L_best] [n_iter] [n_steps_between_global_updates]")
    exit(1)

# Assign input parameters
tau_0 = float(sys.argv[1])
rho = float(sys.argv[2])
L_best = float(sys.argv[3])
n_iter = int(sys.argv[4])
n_steps_between_global_updates = int(sys.argv[5])

print("Using tau_0 = %f, rho = %f, L_best = %f, n_iter = %d, n_steps_between_global_updates = %d" % (tau_0, rho, L_best, n_iter, n_steps_between_global_updates))

# Initialize pheromone for step zero
tau_vs_iteration = np.zeros(n_iter)
tau_vs_iteration[0] = tau_0
iterations = np.arange(0, n_iter)

# Loop over updates
for i_iter in iterations[1:]:

    # Do local update on every round
    tau_vs_iteration[i_iter] = (1 - rho) * tau_vs_iteration[i_iter-1] + rho*tau_0
    
    # Do global update as dictated by the user
    if i_iter % n_steps_between_global_updates == 0:    
        tau_vs_iteration[i_iter] = (1 - rho) * tau_vs_iteration[i_iter] + rho / L_best


# Plot the results

plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
rc('text', usetex=True)
plt.rc('font', size = parameters.default_font_size)
plt.plot(iterations, tau_vs_iteration, 'r-', linewidth=3)
plt.plot(iterations, np.ones(n_iter)*tau_0)
plt.xlabel('Iteration')
plt.ylabel('Pheromone concentration')
plt.legend(['$\\tau$', '$\\tau_0$'])
plt.show()

exit(0)
