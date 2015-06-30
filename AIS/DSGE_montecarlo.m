outfile = "ThetahatsLocalLinear50.out";
mc_reps = 500;

% first design used 10 iters, nparticles = 600
iters = 10;
mixture = 0; % proportion sampled from original prior
initialparticles = 300; % number to take from sample from prior
nparticles = 300; % number per round
particlequantile = 20; % keep the top % of particles
nparticles2 = 600; % number for last round
verbose = false;
AISdraws = 5000; # number of draws from final AIS density
nneighbors = 300;

% design
parameters; % loaded from Common to ensure sync with Gendata

lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
prior_params = [lb ub];
theta0 = lb_param_ub(:,2); % original form
nparams = rows(theta0);

% which statistics to use
load selected; % selected statistics
asbil_selected = selected;

setupmpi; % sets comm world, nodes, node, etc.
asbil_theta = theta0; setupdynare; % sets structures and RNG for simulations
MPI_Barrier(CW);
warning ( "off") ;

% number of particles for each node
particles_per_node = floor(nparticles/(nodes-1));
particles_per_node2 = floor(nparticles2/(nodes-1));

%  frontend: load data and make containers
if !node
	% here, you need to provide code that defines USERthetaZ,
	% which is a reasonably large number set of
	% [theta  Z] where theta is a draw from prior
	% and Z is the output of aux_stat
	load simdata.paramspace;
	USERthetaZ = clean_data(simdata);
	% containers
	thetahats = zeros(mc_reps, nparams);
	stds = zeros(mc_reps, nparams);
endif

for rep = 1:mc_reps
	% the 'true' Zn
	if node==1 % simulate on node1 (can't do it on 0, no *.modfile
		asbil_theta = theta0;
		ok = false;
		while !ok
			USERsimulation;
			Zn = aux_stat(data);
			ok = Zn(:,1) != -1000;
		endwhile	
		Zn = Zn(asbil_selected,:);
		MPI_Send(Zn, 0, mytag, CW);
		for i = 2:nodes-1
			MPI_Send(Zn, i, mytag, CW);
		endfor	
	else % receive it on the other nodes
		Zn = MPI_Recv(1, mytag, CW);
	endif	
	Zn = Zn';
	
	% call the algoritm that gets AIS particles
	AIS; # this gets the particles

    % now draw from the AIS density
    reps_per_node = round(AISdraws/(nodes-1));
    AIS2; # this samples from AIS density

	% see the results
	if !node
		thetas = contribs(:,1:nparams);
        Zs = contribs(:, nparams+1:end);
		test = sum(Zs,2) != 0;
        thetas = thetas(test,:);
        Zs = Zs(test,:);
        Z = [Zn; Zs];
        Z2 = Z;
        q = quantile(Z2,0.99);
	    test = Z < q;
        Z2 = test.*Z2 + (1-test).*q;
	    q = quantile(-Z2,0.99);
	    test = -Z2 < q;
	    Z2 = test.*Z2 - (1-test).*q;
        stdZ = std(Z2);
        Z = Z./stdZ;
		Zs = Z(2:end,:);
		Zn = Z(1,:);

        % find neighbors
        [nn_idx, dd] = nearest_neighbors(Zn, Zs, nneighbors);
        % particles with positive weight
        thetas = thetas(nn_idx,:);
        Zs = Zs(nn_idx,:);
        % weights
        AISweights = prior(thetas) ./ AIS_density(thetas, particles(:,1:nparams));
        weight = dd';
        if max(weight) > 0
            weight = 2*weight/max(weight);
        else
            weight = ones(size(weight));
        endif
        weight = AISweights.*normpdf(weight); # AIS_weights != 1 is for SBIL by AIS
        weight = weight/sum(weight(:));
        % local linear regression
        X = [ones(rows(Zs),1) Zs];
        XX = diag(weight)*X;
        b = inv(X'*XX)*XX'*thetas;
        thetahat = [1 Zn]*b;
		% store results
        thetahats(rep,:) = thetahat;
		FN = fopen (outfile, "a");
		fprintf(FN, "%f ", thetahat);
		fprintf(FN, "\n");
		fclose(FN);
		system('sync');
		if rep > 1
			contrib = thetahats(1:rep,:);
			m = mean(contrib);
			s = std(contrib);
			e = contrib - repmat(theta0',rows(contrib),1);
			b = mean(e);
			e = e.^2;
			mse = mean(e);
			rmse = sqrt(mse);
			lb = lb_param_ub(:,1);
			ub = lb_param_ub(:,3);
			priormean = (ub+lb)'/2;
			priorsdev = sqrt(((ub-lb).^2)/12);
			priorsdev = priorsdev';
			priorbias = priormean - theta0';
			priorrmse = sqrt(priorbias.^2 + priorsdev.^2);
			mae = mean(abs(e));
			clabels = char("true", "mean", "pmean", "sdev.","psdev", "bias", "pbias","rmse", "prmse");
			rlabels = char(
			"alpha",
			"beta",
			"delta",
			"gam",
			"rho1",
			"sigma1",
			"rho2",
			"sigma2",
			"nss"
			);
			printf("\n\nSBIL estimation results: rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
		endif
	endif
endfor

if not(MPI_Finalized) MPI_Finalize; endif

