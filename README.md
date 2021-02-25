# HierarchicalMatrixCompletion

Code for the Hierarchical Matrix completion for Inverse Problems

This is a prototype code for performing hierarchical Matrix completion in the context of inverse problems. 

This project has several componenets:

- FEM for the EIT problem using P1 finite elements
- Routines to compute the DtN map, for a given conductivity
- Computation of the derivatives for a typical L^2 misfit function
- Full optimization loop (using fminunc) for solving the EIT
- Code for performing hierarchical matrix completion, using CVX as a backend.

For this code you will need the following dependencies:
- distmesh
- cvx

