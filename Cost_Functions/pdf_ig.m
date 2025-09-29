function [pdf_out] = pdf_ig(x,mu,lam)
% Calculates PDF of first accumulation period

sq = sqrt(lam./(2.*pi.*x.^3));
exp1 = exp(-(lam*(x-mu).^2)./(2.*x.*mu.^2));
pdf_out = sq.*exp1;
end