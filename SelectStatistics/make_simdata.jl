# this generates the parameter draws and
# candidate statistics for the simple
# linear regression example
n = 100
reps = 20500 # 10000 in sample, 500 out of sample
simdata = zeros(reps, 35+6)
for rep = 1:reps
    # draw the regressors
    x = randn(n,4)
    z = [ones(n,1) x]
    # draw the parameters from prior
    b = 4.*rand(5,1) - 2.0
    sig = 5.*rand()
    # generate dependent variable
    y = z*b + sig*randn(n,1)
    # linear model
    bhat1 = z\y
    uhat = y-z*bhat1
    sighat1 = sqrt(uhat'*uhat/(n-size(z,2)))
    # quadratic model
    z = [ones(n,1) x 0.1.*x.^2.]
    bhat2 = z\y
    uhat = y-z*bhat2
    sighat2 = sqrt(uhat'*uhat/(n-size(z,2)))
    # cubic model
    z = [ones(n,1) x 0.1.*x.^2. 0.01.*x.^3.]
    bhat3 = z\y
    uhat = y-z*bhat3
    sighat3 = sqrt(uhat'*uhat/(n-size(z,2)))
    # pure noise
    z = randn(1,5)
    # assemble
    simdata[rep,:] = [b' sig bhat1' sighat1 bhat2' sighat2 bhat3' sighat3 z]
end

writedlm("simdata100", simdata)
