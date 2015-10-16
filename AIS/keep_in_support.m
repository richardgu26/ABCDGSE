function theta = keep_in_support(theta)
		s = size(theta);
        theta = theta(:);
        parameters;
        lb = lb_param_ub(:,1);
        ub = lb_param_ub(:,3);
        theta = max(theta,lb);
        theta = min(theta,ub);
        theta = reshape(theta,s);
endfunction


