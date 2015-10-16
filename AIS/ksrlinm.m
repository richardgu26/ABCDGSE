function r=ksrlin(x,y,x0,weights, confidencelevel)
% KSRLIN   Local linear kernel smoothing regression, extended to vector x
% x, x0: S by d
% y: S by 1
% h: scalar
r.S=length(y);
r.h=h;

d = size(x,2);
k = size(y,2);

diff = x - x0;
z = diff;

r.weights = weights;

r.localconstantmean = y'*(r.weights)/sum(r.weights);

r.localconstantmean = r.localconstantmean';

diff1 = [ones(size(diff,1),1) diff];
y1 = y.*(sqrt(r.weights)*ones(1,k));
x1 = diff1.*(sqrt(r.weights)*ones(1, d+1));
r.locallinearmean = lscov(x1, y1);
y2 = y.*(r.weights*ones(1,k));
x2 = diff1.*(r.weights*ones(1, d+1));
x3 = x2(:,1);

r.locallinearmedian = zeros(d+1,k);
r.localconstantmedian = zeros(1,k);

for pardim = 1:k;
y2_sub = y2(:,pardim);
r.locallinearmedian(:,pardim) = rq(x2, y2_sub, 0.5); 
end;

for pardim = 1:k;
y2_sub = y2(:,pardim);
r.localconstantmedian(:,pardim) = rq(x3, y2_sub, 0.5); 
end;

tau1 = (1-confidencelevel)/2;
tau2 = confidencelevel + tau1;

r.locallinearlowerquantile = zeros(d+1,k); 
for pardim = 1:k;
y2_sub = y2(:,pardim);
r.locallinearlowerquantile(:,pardim) = rq(x2, y2_sub, tau1); 
end;

r.locallinearupperquantile = zeros(d+1,k); 
for pardim = 1:k;
y2_sub = y2(:,pardim);
r.locallinearupperquantile(:,pardim) =rq(x2, y2_sub, tau2); 
end;

r.localconstantlowerquantile = zeros(1,k); 
for pardim = 1:k;
y2_sub = y2(:,pardim);
r.localconstantlowerquantile(:,pardim)=rq(x3, y2_sub, tau1); 
end;

r.localconstantupperquantile = zeros(1,k); 
for pardim = 1:k;
y2_sub = y2(:,pardim);
r.localconstantupperquantile(:,pardim)=rq(x3, y2_sub, tau2); 
end;

r.locallinearmean = r.locallinearmean(1,:);
r.locallinearmedian = r.locallinearmedian(1,:);
r.locallinearlowerquantile = r.locallinearlowerquantile(1,:);
r.locallinearupperquantile = r.locallinearupperquantile(1,:);


