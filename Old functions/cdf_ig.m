function [cdf_out] = cdf_ig(x,mu,lam)
sq = sqrt(lam./x);
xm = x./mu;
phi1 = normcdf(sq.*(xm-1));
exp1 = exp(2.*lam./mu);
phi2 = normcdf((-sq).*(xm+1));

cdf_out = phi1+exp1.*phi2;

% % sq = sqrt(lam/x);
% % denom = sqrt(2)*mu;
% % erfc1 = erfc(sq*(-x+mu)./denom);
% % erfc2 = erfc(sq.*(x+mu)./denom);
% % exp1 = exp(2*lam./mu);
% % 
% % cdf_out = 0.5*erfc1+0.5*exp1*erfc2;

% sq = sqrt(lam./x);
% xm = sqrt(2)*x/mu;
% phi1 = normcdf(x,-sq.*(xm+1));
% exp1 = exp(2*lam./mu);
% phi2 = normcdf(x,(-sq).*(xm-1));
% 
% cdf_out = phi1+exp1.*phi2;


end