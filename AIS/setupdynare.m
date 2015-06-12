UseDynare;
#if node
        % break into pieces
        alpha = asbil_theta(1,:);
        beta  = asbil_theta(2,:);
        delta = asbil_theta(3,:);
        gam   = asbil_theta(4,:);
        rho1   = asbil_theta(5,:);
        sigma1 = asbil_theta(6,:);
        rho2   = asbil_theta(7,:);
        sigma2 = asbil_theta(8,:);
        nss   = asbil_theta(9,:);
        save parameterfile  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;

        % get the RNG state on this node, to re-establish separate
        % states on nodes after running Dynare, which synchronizes them
        RNGstate = rand('state');

        % solve model once on each node, to get Dynare structures ready
        % for simulations
        command = sprintf("dynare SimpleModel%d noclearall", node);
        eval(command);

        % re-set the seed
        ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
        set_dynare_seed(ss);
#end
