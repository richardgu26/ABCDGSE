#= 
MPI version of selection code for the simple example in
Creel and Kristensen (2015) "On Selection of Statistics for Approximate
Bayesian Computing (or the Method of Simulated Moments)"

to run this, enter mpirun -np X julia Select_mpi.jl
from the system prompt. X should be an even divisor of 
the number of reps, and should be less than the total
number of available CPU cores
=#

using Distances
using Distributions
import MPI
include("SelectionAlgorithm.jl")
function main()

blas_set_num_threads(2)
if ~MPI.Initialized() MPI.Init() end
comm = MPI.COMM_WORLD
MPI.Barrier(comm)
node = MPI.Comm_rank(comm)
nodes = MPI.Comm_size(comm)

    # run with -np X+1, X should be an even divisor of reps
    reps = 25 # how many repetitions of algorithm

    nparams = 9 # number of parameters in data set
    S = 30000  # size of paramspace sample
    S2 = 1000  # size of test sample
 
    # data file to read, initial nparams columns are parameters,
    # rest are candidate statistics
    infile = "simdata.dsge"
    # output file
    outfile = "battery.dsge"
    whichdep = 0; # set to 0 for all parameters, or target to a specific parameter

    # read the data, and split into params and stats
    simdata = readdlm(infile) 
    theta = simdata[1:S+S2,1:nparams]
    Z = simdata[1:S+S2,nparams+1:end]
    simdata = 0. # save some memory

    # standardize and normalize stats and parameters
    # Z is divided by trimmed std. dev, so that outliers
    # don't compress normal draws
    Z2 = copy(Z)
    dimZ = size(Z,2)
    @inbounds for i = 1:dimZ
        q = quantile(Z2[:,i],0.99)
        # top bound
        test =  Z2[:,i] .< q
        Z2[:,i] = Z2[:,i].*test  + q.*(1. - test)
        q = -quantile(-Z2[:,i],0.99)
        # bottom bound
        test =  Z2[:,i] .> q
        Z2[:,i] = Z2[:,i] .* test + q.*(1. - test)
    end
    stdZ = std(Z2,1) 
    Z = (Z .-mean(Z))./stdZ 
    theta = theta .- mean(theta)
    s = std(theta,2)
    s = s.*(s.>0.) .+ (s.<=0) # allows for fixed parameters
    theta = theta ./s
    # split into in and out of sample
    Z_in = Z[1:S,:]
    Z_out = Z[S+1:S+S2,:]
    theta_in = theta[1:S,:]
    theta_out = theta[S+1:S+S2,:]
 
    # number of neighbors for SBIL fits
    n, dimZ = size(Z_in)
    #neighbors = round(Int, floor(1.*n^0.25))
    neighbors = 100  # local linear (now the default) needs more neighbors
    # SA controls
    temperature = 0.005 # 0.02 works well, given that parameters are standardized and normalized
    r = 1. # factor to set nt
    nt = r*dimZ # number of evaluations between temp reductions
    rt = 0.75 # rate of temp reduction
    r = 10; # total evals is this number times dimZ
    maxevals = r*dimZ
    # do sa
    if node>0
        for i = 1:round(Int,reps/(nodes-1))
            selected, obj_value = sa_for_selection(Z_in, Z_out, theta_in, theta_out, temperature, nt, rt, maxevals, neighbors, whichdep)
	        MPI.Isend([obj_value vec(selected)'], 0, node, comm)
        end    
    else # frontend
        done = false
        received = 0 # this many ranks need to report
        results = zeros(reps, dimZ+1)
        while received < reps
            # retrieve results from compute nodes
            sleep(0.01) # this takes a while, may as well not flood
            for i = 1:nodes-1
                # compute nodes have results yet?
                ready = false
                ready, junk = MPI.Iprobe(i, i, comm) # check if message pending
                if ready
                    # get it if it's there
                    contrib = zeros(1,dimZ+1)
                    MPI.Recv!(contrib, i, i, comm)
                    received +=1
                    results[received,:] = contrib
                end
            end
        end
        writedlm(outfile, results) 
    end 
    # clean up MPI communicator
    MPI.Barrier(comm)
    if ~(MPI.Finalized()) MPI.Finalize() end
end
# call the main function
main()


