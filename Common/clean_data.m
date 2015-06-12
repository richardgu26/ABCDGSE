function d = clean_data(simdata)

%	printf("%d rows read in\n", rows(simdata));
	test = any(isnan(simdata'));
	test = (test ==0);
	simdata = simdata(test,:);
	test = any(isinf(simdata'));
	test = (test ==0);
	simdata = simdata(test,:);
%	printf("%d rows without NaN and Inf\n", rows(simdata));
	test = simdata == -1000;
	test = sum(test,2);
	test = (test ==0);
	simdata = simdata(test,:);
%	printf("%d rows after bad data code \n", rows(simdata));

	d = simdata;
end
