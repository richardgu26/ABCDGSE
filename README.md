# ABCDGSE
This directory holds the files for replicating the procedure used
to get results for the DSGE model in the paper "Indirect Likelihood
Inference" by Michael Creel and Dennis Kristensen

The estimator is an Approximate Bayesian Computing estimator,
computed using importance sampling and local linear nonparametric
regression. 

Requires Open MPI, Julia (with MPI and Distances packages), Octave (with MPI
package) and other supporting code available at https://github.com/mcreel/Econometrics

To replicate the results, execute "sh MasterScript" from the bash prompt.

For questions, write to michael.creel@uab.es
