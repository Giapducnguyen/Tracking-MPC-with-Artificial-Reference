# Tracking-MPC-with-Artificial-Reference
This repository contains the file reproducing the simulation results of the paper: Limón, Daniel, et al. "MPC for tracking piecewise constant references for constrained linear systems." Automatica 44.9 (2008): 2382-2387.

Perhaps the most confusing part is how to compute the set $\mathcal{O}_{\infty, \lambda}^{w}$. For that reason, it is recommended to read the following paper:

Gilbert, Elmer G., and K. Tin Tan. "Linear systems with state and control constraints: The theory and application of maximal output admissible sets." IEEE Transactions on Automatic control 36.9 (1991): 1008-1020.

It is noted that MATLAB function $\texttt{fmincon}$ was used for fast prototyping. However, QP solvers can be used instead since the problem is linear. Also, set operations were performed using the MPT3 toolbox available at https://www.mpt3.org/  
