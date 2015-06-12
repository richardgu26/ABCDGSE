## Copyright (C) 2014 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author: Michael Creel <michael@pelican>
## Created: 2014-10-22

%% this is the main asbil algorith code, now in one place for synchronization

% first define functions, then the main algorithm
1;

function theta_s = sample_from_AIS(particles)
    delta = 0.1*std(particles);    
	i = randi(rows(particles));
	theta_s = particles(i,:);
	theta_s = theta_s + delta .*randn(size(delta));
	theta_s = theta_s';
endfunction


# the importance sampling density: mixture of normals
function dens = AIS_density(thetas, particles)
    # what scaling to use here?
    delta  = 0.1*std(particles,1);
    delta = delta.^2; # variances
    sig = diag(delta);
    nparticles = size(particles,1);
    dens = zeros(rows(thetas),1);
    for i = 1:nparticles
        mu = particles(i,:);
        dens += mvnpdf(thetas, mu, sig);
    end
    dens = dens/nparticles;
end    

% main algorithm code, kept in one place to ensure Monte Carlo and estimation are doing the same thing
if node
    thetas = zeros(reps_per_node, nparams);
    Zs = zeros(reps_per_node, size(Zn,2));
	for i = 1:reps_per_node;	
        thetass = sample_from_AIS(particles(:,1:nparams));
        thetas(i,:) = thetass';
		ok = all((thetass >= lb) & (thetass <= ub),1);
        if ok # if support conditions ok, do the simulation, otherwise, no
            ok = false;
            while !ok # repeat until good draw obtained
                asbil_theta = thetass;
                USERsimulation;
                Z = aux_stat(data);
                ok = Z(:,1) != -1000;
            endwhile	
            Z = Z(asbil_selected,:);
            Zs(i,:) = Z';
        endif
	endfor
    contribs = [thetas Zs];
	MPI_Send(contribs, 0, mytag, CW); % send selected to frontend, to sync across nodes

else % frontend
  		% receive selected particles from nodes
		for i = 1:nodes-1
			contrib = MPI_Recv(i, mytag, CW);
			if (i == 1)
				contribs = contrib;	
			else
				contribs = [contribs; contrib];
			endif
		endfor
endif

# get rid of draws that were out of support
contribs = clean_data(contribs);
