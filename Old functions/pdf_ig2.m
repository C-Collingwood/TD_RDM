function [pdf_out] = pdf_ig2(x,b,a,mu1,s)
mu = @(b)(b./mu1);
lam = @(b)(b/s).^2;

sq = @(b)(sqrt(lam(b)./(2.*pi.*x.^3)));
exp1 = @(b)exp(-(lam(b)*(x-mu(b)).^2)./(2.*x.*mu(b).^2));
pdf_out = integral(@(b)(sq(b).*exp1(b)),b-a,b);


end

