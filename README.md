# ABCDGSE
This directory holds the files for replicating the procedure used
to get results for the DSGE model in the paper "Bayesian Indirect
Inference and the ABC of GMM" by Creel, Gao, Hong and Kristensen, which
will be available shortly. The file DSGE.pdf is an extract, which
gives a description and the results.

The estimator is an Approximate Bayesian Computing estimator,
computed using importance sampling and local linear nonparametric
regression. 

Requires Open MPI, Julia (with MPI and Distances packages), Octave (with MPI
package) and other supporting code available at https://github.com/mcreel/Econometrics

To replicate the results, execute "sh MasterScript" from the bash prompt.

For questions, write to michael.creel@uab.es
