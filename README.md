# kmedian

Illustration of Lagrangian techniques, as applied to the k-median problem. 

Includes code for:
1. Creating random test instances
2. Solving k-median IP (and LP relaxation) via Gurobi
3. Solving k-median with branch-and-bound (using Lagrangian bounds)

The Lagrangian outer problem is solved with Eugene Lykhovyd's implementation of Shor's r-algorithm (version 21-May-2019)
  https://github.com/zhelih/ralg
  
Version 21-May-2019 of Lykhovyd's code requires an MKL installation:
  https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-2018-install-guide
  
See Lykhovyd's makefile for the necessary includes/libraries. The 5-June-2019 version of Lykhovyd's code does not require MKL and will be slower when MKL is not used.
