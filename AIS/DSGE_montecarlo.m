outfile = "junk.out";
mc_reps = 1000; % number of MC reps
nworkers = 20;  % number of worker MPI ranks

% controls for creating the adaptive importance sampling density
iters = 10;
initialparticles = nworkers*round(300/nworkers); % number to take from sample from prior
nparticles = nworkers*round(300/nworkers); % number per round
particlequantile = 20; % keep the top % of particles
verbose = false;

% controls for drawing the final sample from mixture of AIS and prior
mixture = 0.5; % proportion sampled from original prior 
AISdraws = nworkers*round(5000/nworkers); # number of draws from final AIS density

% controls for the nonparametric fits
nneighbors = 500;    

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

%  frontend: load data and make containers
if !node
	% here, you need to provide code that defines USERthetaZ,
	% which is a reasonably large number set of
	% [theta  Z] where theta is a draw from prior
	% and Z is the output of aux_stat
	load simdata.paramspace;
	USERthetaZ = clean_data(simdata);
	% containers
	thetahatsLC = zeros(mc_reps, nparams);
	thetahatsLC50 = zeros(mc_reps, nparams);
	thetahatsLL = zeros(mc_reps, nparams);
	thetahatsLL50 = zeros(mc_reps, nparams);
	thetahatsLQ = zeros(mc_reps, nparams);
	thetahatsLQ50 = zeros(mc_reps, nparams);
	cilower = zeros(mc_reps, nparams);
	ciupper = zeros(mc_reps, nparams);
endif

for rep = 1:mc_reps
    % the 'true' Zn
    if node==1 % simulate on node1 (can't do it on 0, no *.modfile
        asbil_theta = theta0;
        ok = false;
        while !ok    
            USERsimulation;
            Zn = aux_stat(data);
            ok = Zn(1,:) != -1000;
        endwhile	
        Zn = Zn(asbil_selected,:);
        for i = 2:nodes-1
            MPI_Send(Zn, i, mytag, CW);
        endfor	
        MPI_Send(Zn, 0, mytag, CW);
    else % receive it on the other nodes
        Zn = MPI_Recv(1, mytag, CW);
        if  !node
        endif    
    endif
    MPI_Barrier(CW);    
    Zn = Zn';
 
	% call the algoritm that gets AIS particles
	if !node
            printf("starting AIS\n");
            tic;
    endif  
    CreateAIS; # this gets the particles

    % now draw from the AIS density
    reps_per_node = round(AISdraws/(nodes-1));
   	if !node
            toc;
            printf("starting AIS2\n");
            tic;
    endif  
    SampleFromAIS; # this samples from AIS density


	% see the results

   	if !node
        toc;    
        % create the bandwidths used in tuning
        nbw = 30;
        bandwidths = zeros(nbw,1);
        for i = 1:nbw
                bandwidths(i,:) = 0.1 + (10-0.1)*((i-1)/(nbw-1))^2;
        endfor      
        % selected bws from tuning
        % selected using prior
        %bwselect = [6 6 6 7 7 7 7 8 6];
        %bwselectCI = [ 12 22 14 12 15 12 15 11 21 ];
        
        % selected using local
        bwselect = [6 4 5 8 6 6 7 7 5];
        bwselectCI = [13 27 17 3 24 12 11 10 12];
        
        bandwidthsCI = bandwidths(bwselectCI,:);
        bandwidths = bandwidths(bwselect,:);

	    printf("starting fit and CI\n");
        tic;
		thetas = contribs(:,1:nparams);
        Zs = contribs(:, nparams+1:end);
        test = sum(Zs,2) != 0;
        thetas = thetas(test,:);
        Zs = Zs(test,:);
        Z = [Zn; Zs];
        
        % first pre-whiten using all draws
        %q = quantile(Z,0.99);
        %test = Z < q;
        %Z = test.*Z + (1-test).*q;
        %q = quantile(-Z,0.99);
        %test = -Z < q;
        %Z = test.*Z - (1-test).*q;
	    stdZ = std(Z);
        Z = Z ./stdZ;
        Zs = Z(2:end,:);
		Zn = Z(1,:);

        %AISweights = prior(thetass) ./ (mixture*prior(thetass) +(1-mixture)*AIS_density(thetass, particles(:,1:nparams)));
        AISweights = 1;
        weights = zeros(rows(Zs),9);
        for i = 1:9
            weights(:,i) = __kernel_normal((Zs-Zn)/bandwidths(i,:));
        endfor    
        weights = AISweights.*weights; # AIS_weights != 1 is for SBIL by AIS
        weights = weights./sum(weights);
        thetahatLC = zeros(1,9);
        thetahatLL = zeros(1,9);
        thetahatLQ = zeros(1,9);
        thetahatLC50 = zeros(1,9);
        thetahatLL50 = zeros(1,9);
        thetahatLQ50 = zeros(1,9);
        for i = 1:9
                r = LocalConstant(thetas(:,i), weights(:,i), false);
                thetahatLC(:,i) = r.mean;
                thetahatLC50(:,i) = r.median;
                r = LocalPolynomial(thetas(:,i), Zs, Zn, weights(:,i), false);
                thetahatLL(:,i) = r.mean;
                thetahatLL50(:,i) = r.median;
                r = LocalPolynomial(thetas(:,i), Zs, Zn, weights(:,i), false, 2);
                thetahatLQ(:,i) = r.mean;
                thetahatLQ50(:,i) = r.median;
        endfor
        thetahatLC = keep_in_support(thetahatLC);
        thetahatLC50 = keep_in_support(thetahatLC50);
        thetahatLL = keep_in_support(thetahatLL);
        thetahatLL50 = keep_in_support(thetahatLL50);
        thetahatLQ = keep_in_support(thetahatLQ);
        thetahatLQ50 = keep_in_support(thetahatLQ50);
        % confidence intervals
        % weights
        %AISweights = prior(thetas) ./ (mixture*prior(thetas) + (1-mixture)*AIS_density(thetas, particles(:,1:nparams)));
        AISweights = 1;
        weights = zeros(rows(Zs),9);
        for i = 1:9
            weights(:,i) = __kernel_normal((Zs-Zn)/bandwidthsCI(i,:));
        endfor         
        weights = AISweights.*weights; # AIS_weights != 1 is for SBIL by AIS
        %weights = weights + eps; # keep positive
        weights = weights./sum(weights);
        % for CIs, use local constant
        lower = zeros(9,1);
        upper = lower;
        for i = 1:9
            r = LocalConstant(thetas(:,i), weights(:,i), true);
            lower(i,:) = r.c;
            upper(i,:) = r.d;
        endfor    
		% results
        if rep == 1
                in_ci = zeros(rows(theta0),1);
        endif
        % CI coverage
        in10 = ((theta0 > lower) & (theta0 < upper));
        in_ci = in_ci + in10;
        thetahatsLC(rep,:) = thetahatLC;
        thetahatsLC50(rep,:) = thetahatLC50;
        thetahatsLL(rep,:) = thetahatLL;
        thetahatsLL50(rep,:) = thetahatLL50;
        thetahatsLQ(rep,:) = thetahatLQ;
        thetahatsLQ50(rep,:) = thetahatLQ50;
        cilower(rep,:) = lower';
        ciupper(rep,:) = upper';
        save(outfile, "thetahatsLC", "thetahatsLC50","thetahatsLL", "thetahatsLL50", "thetahatsLQ", "thetahatsLQ50", "cilower", "ciupper");
		if rep > 1
			contrib = thetahatsLC(1:rep,:);
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
			printf("\n\nEstimation results (LC mean): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
            % now median
			contrib = thetahatsLC50(1:rep,:);
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
			printf("\n\nEstimation results (LC median): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
			contrib = thetahatsLL(1:rep,:);
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
			printf("\n\nEstimation results (LL mean): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
            % now median
			contrib = thetahatsLL50(1:rep,:);
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
			printf("\n\nEstimation results (LL median): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
			contrib = thetahatsLQ(1:rep,:);
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
			printf("\n\nEstimation results (LQ mean): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
            % now median
			contrib = thetahatsLQ50(1:rep,:);
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
			printf("\n\nEstimation results (LQ median): rep %d\n", rep);
			prettyprint([theta0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
			printf("\n\n");
            printf("90%% CI coverage: \n");
            disp(in_ci/rep);
            disp([lower upper]);
            printf("\n");
		endif
    endif
endfor

if not(MPI_Finalized) MPI_Finalize; endif

