function [cdf_out] = cdf_ig2(x,b,a,mu1,s)
mu = @(b)(b./mu1);
lam = @(b)(b/s).^2;
sq = @(b)sqrt(lam(b)./x);
xm = @(b)(x./mu(b));
phi1 = @(b)(sq(b).*(xm(b)-1));
exp1 = @(b)exp(2.*lam(b)./mu(b));
phi2 = @(b)normcdf((-sq(b)).*(xm(b)+1));




cdf_out = integral(@(b)(phi1(b)+exp1(b).*phi2(b)), b-a,b);

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