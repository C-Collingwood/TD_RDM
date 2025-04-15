function [cdf_out] = cdf_ig(x,mu,lam)
sq = sqrt(lam./x);
xm = x./mu;
phi1 = normcdf(sq.*(xm-1));
exp1 = exp(2.*lam./mu);
phi2 = normcdf((-sq).*(xm+1));

cdf_out = phi1+exp1.*phi2;

end