# kmedian

Illustration of Lagrangian techniques, as applied to the k-median problem. 

Includes code for:
1. Creating random test instances
2. Solving k-median IP (and LP relaxation) via Gurobi
3. Solving k-median with branch-and-bound (using Lagrangian bounds)

The Lagrangian outer problem is solved with Eugene Lykhovyd's implementation of Shor's r-algorithm
  https://github.com/zhelih/ralg
