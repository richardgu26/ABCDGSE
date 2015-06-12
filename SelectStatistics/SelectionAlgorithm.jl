#=
These are the functions for selection of statistics:
* knn_regression: simple KNN smoother
* select_obj: the objective function to be minimized.
* sa_for_selection: the SA algorithm, adapted for subset selection
=#

# this is a simple uniform kernel KNN smoother
function knn_regression(y, x, xeval, neighbors)
    n = size(xeval,1)
    nn, k = size(y)
    yhat = zeros(n,k)
    ind = zeros(Int64,nn)
    distances = pairwise(Euclidean(),x', xeval') # get all distances
    for i = 1:n
        ind = sortperm(sub(distances,:,i)) # indices of k nearest neighbors
        @inbounds yhat[i,:] = mean(y[ind[1:neighbors],:],1) # fit is mean of nearest neighbors
    end
    return yhat
end    


# local linear KNN smoother
function knn_ll_regression(y, x, xeval, neighbors)
    n = size(xeval,1)
    nn, k = size(y)
    yhat = zeros(n,k)
    ind = zeros(Int64,nn)
    distances = pairwise(Euclidean(),x', xeval') # get all distances
    @inbounds for i = 1:n
        di = sub(distances,:,i)
        ind = sortperm(di) # indices of k nearest neighbors
        selected = vec(ind[1:neighbors,:])
        yi = y[selected,:]
        xi = x[selected,:]
        di2 = di[selected,:]
        m = maximum(di2)
        if m > 0
            weight = 2.*di2/m
        else
            weight = 0.
        end    
        weight = pdf(Normal(),weight)
        weight = weight/sum(weight)
        X = [ones(neighbors,1) xi]
        XX = weight .* X;
        b = inv(X'*XX)*XX'*yi
        yhat[i,:] = [1. xeval[i,:]]*b
    end
    return yhat
end    

# the objective function to be minimized: MAE, plus a penalty
function select_obj(selected, Z_in, Z_out, theta_in, theta_out, k, whichdep)
    s2 = vec(map(Bool,selected)) # SA alg. works with numbers, here we want boolean
    Z_in1 = Z_in[:,s2]
    Z_out1 = Z_out[:,s2]
	if sum(s2) != 0
		thetahat = knn_ll_regression(theta_in, Z_in1, Z_out1, k)
		if whichdep == 0  # ordinary, for all parameters
            obj_value = mean(abs(theta_out-thetahat))
        else # targeted to certain parameter
            mae = mean(abs(theta_out[:,whichdep]-thetahat[:,whichdep]))  # 3 is delta
		    obj_value =  mae*(1. + 0.05*sum(s2)) # objective is MAE, plus 5% penalty per statistic
        end
    else
		obj_value = 1000.
	end
    return obj_value
end

# the simulated annealing algorithm
function sa_for_selection(Z_in, Z_out, theta_in, theta_out, temperature, nt, rt, maxevals, k, whichdep=0)
    # dimension of statistics
   	dimZ = size(Z_in,2)
   	# Initial trial value and obj. function 
    x = ones(dimZ,1)
    f = select_obj(x, Z_in, Z_out, theta_in, theta_out, k, whichdep)
    xopt = copy(x)
	fopt = copy(f) # give it something to compare to
	# main loop
    for fevals = 1:maxevals
		h = rand(1:dimZ) # which statistic to add or remove
		xp = copy(x)
		@inbounds xp[h,:] = 1.0-x[h,:] # switch from 0->1 or 1->0
		# Evaluate function at new point
        fp = select_obj(xp, Z_in, Z_out, theta_in, theta_out, k, whichdep)
		fevals += 1
		#  Accept the new point if the function value decreases
		if (fp <= f)
            x = copy(xp)
			f = copy(fp)
			#  If lower than any other point, record as new optimum
			if(fp < fopt)
				xopt = copy(xp)
                x = copy(xp)
				fopt = copy(fp)
            end
		# If the point is higher, use the Metropolis criterion
		else
			p = exp(-(fp - f) / temperature)
			if (rand() < p)
		        x = copy(xp)
				f = copy(fp)
            end
        end
        # reduce temperature every nt evaluations
        # and set current point to best so far
	    if mod(fevals,nt)==0.
            temperature = rt * temperature
            x = copy(xopt)
            f = copy(fopt)
        end    
        # intermediate results
	    if mod(fevals,1)==0.
            println()
            println("fevals: ", fevals, " temperature: ", temperature)
		    println("best obj so far: ", fopt, " number selected: ", sum(xopt))
            a = 1:size(xopt,1)
            selected = vec(map(Bool,xopt))
            a = a[selected]
            println("selected: ", a)
        end    
     end
    return xopt, fopt
end    


